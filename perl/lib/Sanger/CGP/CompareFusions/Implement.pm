package Sanger::CGP::CompareFusions::Implement;
##########LICENCE ##########
#Copyright (c) 2015 Genome Research Ltd.
###
#Author: Cancer Genome Project <cgpit@sanger.ac.uk>
###
#This file is part of cgpRna.
###
#cgpRna is free software: you can redistribute it and/or modify it under
#the terms of the GNU Affero General Public License as published by the
#Free Software Foundation; either version 3 of the License, or (at your
#option) any later version.
###
#This program is distributed in the hope that it will be useful, but
#WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
#General Public License for more details.
###
#You should have received a copy of the GNU Affero General Public
#License along with this program. If not, see
#<http://www.gnu.org/licenses/>.
###
#1. The usage of a range of years within a copyright statement contained
#within this distribution should be interpreted as being equivalent to a
#list of years including the first and last year specified and all
#consecutive years between them. For example, a copyright statement that
#reads ‘Copyright (c) 2005, 2007- 2009, 2011-2012’ should be interpreted
#as being identical to a statement that reads ‘Copyright (c) 2005, 2007,
#2008, 2009, 2011, 2012’ and a copyright statement that reads ‘Copyright
#(c) 2005-2012’ should be interpreted as being identical to a statement
#that reads ‘Copyright (c) 2005, 2006, 2007, 2008, 2009, 2010, 2011,
#2012’."
##########LICENCE ##########
##########
use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use English qw( -no_match_vars );
use Const::Fast qw(const);
use File::Which qw(which);
use File::Basename;
use FindBin qw($Bin);
use File::Spec;
use PCAP::Cli;
use PCAP::Threaded;
use Sanger::CGP::CompareFusions::FusionAnnotation;
use Sanger::CGP::CgpRna;
use Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource;
use Sanger::CGP::Vagrent::Data::GenomicRegion;
use Sanger::CGP::Vagrent::Data::Transcript;

use Data::Dumper;

const my $BEDTOOLS_CLOSEST => q{ closest -s -a %s -b %s | sort -k4,4 > %s};
const my $BEDTOOLS_PAIRTOPAIR => q{ pairtopair -a %s -b %s -slop 5 > %s};
const my $SORT => q{ sort -k%d,%d %s > %s};
const my $TRI_SORT => q{ sort -k%d,%d -k%d,%dr -k%d,%dr %s > %s};
const my $JOIN => q{ join -1 8 -2 8 %s %s > %s };
const my $OUTPUT_HEADER => "sample\tfusion_name\talgorithm\tconfidence_score\tstar_junction\ttophat_junction\tdefuse_junction\tdefuse_cluster_id\tstar_junction_reads\tstar_spanning_frags\ttophat_junction_reads\ttophat_spanning_frags\tdefuse_splitr_count\tdefuse_span_count\t5'_gene\t5'_gene_id\t5'_chr\t5'_pos\t5'_strand\t3'_gene\t3'_gene_id\t3'_chr\t3'_pos\t3'_strand\t5'_transcript_id\t5'_transcript_src\t5'_exon_num\t5'_exon_start\t5'_exon_end\t3'_transcript_id\t3'_transcript_src\t3'_exon_num\t3'_exon_start\t3'_exon_end\tdefuse_splitr_sequence\ttophat_splitr_sequence\n";

# This filter on biotypes is currently not used in subroutine filter_gtf (uncomment the line to switch on).
my %ALLOWED_BIOTYPES = (
	antisense => 1,
	IG_C_gene => 1,
	IG_D_gene => 1,
	IG_J_gene => 1,
	IG_V_gene => 1,
	lincRNA => 1,
	miRNA => 1,
	processed_transcript => 1,
	retained_intron => 1,
	protein_coding => 1,
	rRNA => 1,
	sense_intronic => 1,
	sense_overlapping => 1,
	snoRNA => 1,
	snRNA => 1,
	TR_C_gene => 1,
	TR_D_gene => 1,
	TR_J_gene => 1,
	TR_V_gene => 1,
);

my %CONFIDENCE_SCORES = (
  STD => '76%',
  SD => '77%',
  ST => '48%',
  TD => '25%',
  D => '58%',
  T => '50%',
  S => '22%',
);

# Position of the columns in the tophat-fusion filtered file used to format the bed file.
const my $TOPHAT_SPLIT_CHAR => '\t';
const my $TOPHAT_GENE1 => 3;
const my $TOPHAT_CHR1 => 4;
const my $TOPHAT_POS1 => 5;
const my $TOPHAT_STRAND1 => 13;
const my $TOPHAT_GENE2 => 6;
const my $TOPHAT_CHR2 => 7;
const my $TOPHAT_POS2 => 8;
const my $TOPHAT_STRAND2 => 14;
const my $TOPHAT_SPAN_READS => 9;
const my $TOPHAT_SPAN_MATE_PAIRS => 10;
const my $TOPHAT_SPAN_MATE_PAIRS2 => 11;
const my $TOPHAT_SCORE => 12;
const my $TOPHAT_BREAKREF => 1;
const my $TOPHAT_HEADER_PATTERN => 'num_spanning_reads';

# Position of the columns in the deFuse output file used to format fusion breakpoint references.
const my $DEFUSE_SPLIT_CHAR => '\t';
const my $DEFUSE_CHR1 => 26;
const my $DEFUSE_POS1 => 39;
const my $DEFUSE_STRAND1 => 36;
const my $DEFUSE_GENENAME1 => 32;
const my $DEFUSE_GENEID1 => 22;
const my $DEFUSE_CHR2 => 27;
const my $DEFUSE_POS2 => 40;
const my $DEFUSE_STRAND2 => 37;
const my $DEFUSE_GENENAME2 => 33;
const my $DEFUSE_GENEID2 => 23;
const my $DEFUSE_BREAKREF => 1;
const my $DEFUSE_CLUSTER_ID => 2;
const my $DEFUSE_SEQUENCE => 3;
const my $DEFUSE_SPLIT_READS => 4;
const my $DEFUSE_SPAN_READS => 62;
const my $DEFUSE_HEADER_PATTERN => 'cluster_id';

# Position of the columns in the star-fusion output file used to format fusion breakpoint references.
const my $STAR_SPLIT_CHAR => '\t';
const my $STAR_FUSION_NAME => 2;
const my $STAR_JUNCTION_READS => 3;
const my $STAR_SPANNING_FRAGS => 4;
const my $STAR_CHR1 => 7;
const my $STAR_POS1 => 8;
const my $STAR_STRAND1 => 9;
const my $STAR_CHR2 => 13;
const my $STAR_POS2 => 14;
const my $STAR_STRAND2 => 15;
const my $STAR_BREAKREF => 1;
const my $STAR_GENENAME1 => 5;
const my $STAR_GENEID1 => 6;
const my $STAR_GENENAME2 => 11;
const my $STAR_GENEID2 => 12;
const my $STAR_DIS_EXON1 => 10;
const my $STAR_DIS_EXON2 => 16;
const my $STAR_HEADER_PATTERN => 'fusion_name';

sub annotate_bed {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $sample = $options->{'sample'};
  my $exon_gtf = filter_gtf($options, 'exon');

  my $break1_file;
  my $break2_file;
	
  opendir(my $dh, $tmp);
  while(my $file = readdir $dh) {
    $break1_file = File::Spec->catfile($tmp, $file) if($file =~ m/^$sample.1.bed/);
    $break2_file = File::Spec->catfile($tmp, $file) if($file =~ m/^$sample.2.bed/);
  }
  closedir($dh);
	
  my $break1_annotated_file = $break1_file;
  my $break2_annotated_file = $break2_file;
  $break1_annotated_file =~ s/bed/ann/;
  $break2_annotated_file =~ s/bed/ann/;
  
  # Format the bedtools closest commands. Use the filtered exon and gene gtf files to ensure we are getting back the relevant features of interest.
  my $prog = _which('bedtools');
  my $command1 = $prog . sprintf $BEDTOOLS_CLOSEST, $break1_file, $exon_gtf, $break1_annotated_file;
  my $command2 = $prog . sprintf $BEDTOOLS_CLOSEST, $break2_file, $exon_gtf, $break2_annotated_file;
	
  my @commands = ($command1,$command2);

  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

  return 1;
}

sub annotation_sort {
        my $a_ccds = 0;
        my $b_ccds = 0;
        if(defined($a->getCCDS) && $a->getCCDS ne ''){
                $a_ccds = 1;
        }
        if(defined($b->getCCDS) && $b->getCCDS ne ''){
                $b_ccds = 1;
        }
        my $ccds_cmp = $b_ccds <=> $a_ccds;
        if($ccds_cmp == 0){
                my $a_cds_len = $a->getCdsLength;
                my $b_cds_len = $b->getCdsLength;
                my $cds_len_cmp = $b_cds_len <=> $a_cds_len;
                if($cds_len_cmp == 0){
                        return $b->getmRNALength <=> $a->getmRNALength;
                } else {
                        return $cds_len_cmp;
                }
        } else {
                return $ccds_cmp;
        }
}

sub check_input {
  my $fusion_file = shift;
	
  # Check the file exists
  PCAP::Cli::file_for_reading('fusion-file', $fusion_file);
	
  # Read the header and identify the file format
  open my $file, '<', $fusion_file; 
  my $firstLine = <$file>;
  close $file;
	
  my $source;
  # Check for tophat format
  if($firstLine =~ m/$TOPHAT_HEADER_PATTERN/){
    $source = "tophat";
  }
  # defuse check
  elsif($firstLine =~ m/$DEFUSE_HEADER_PATTERN/){
    $source = "defuse";
  }
  # star check
  elsif($firstLine =~ m/$STAR_HEADER_PATTERN/){
    $source = "star";
  }
  else{
    die "Unrecognised file type or the file is missing the header record\n";
  }
	
  return $source;
}

sub collate_annotation {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $sample = $options->{'sample'};
  my $annot_file1 = File::Spec->catfile($tmp, "$sample.1.ann_final");
  my $annot_file2 = File::Spec->catfile($tmp, "$sample.2.ann_final");
  my $annot_file_full = File::Spec->catfile($tmp, "$sample.final");
  
  open (my $ifh1, $annot_file1) or die "Could not open file '$annot_file1' $!";
  open(my $ofh1, '>>', $annot_file_full) or die "Could not open file $annot_file_full $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    my @fields = split "\t", $line;
    my @fields2;
    #pop @fields;
    my $breakpoint_id = $fields[0];
    open (my $ifh2, $annot_file2) or die "Could not open file '$annot_file2' $!";
    while (<$ifh2>) {
      chomp;
      my $line2 = $_;
      if($line2 =~ m/^$breakpoint_id\s/){
        @fields2 = split "\t", $line2;
        shift(@fields2);
				shift(@fields2);
				shift(@fields2);
				#push(@fields, @fields2);
				last;
      }
    }
    if(scalar @fields2 > 0){
      pop @fields;
    }
    push(@fields, @fields2);
    close($ifh2);
    print $ofh1 join("\t", @fields)."\n";
  }
  close($ofh1);
  close($ifh1);
  
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  
  return 1;
}

sub create_junction_bedpe {
  my ($index, $options) = @_;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
  my $sample = $options->{'sample'};
  
  my $file = $options->{'fusion_files'}->{$index}->{'name'};
  my $filetype = $options->{'fusion_files'}->{$index}->{'format'};
  my $output1 = File::Spec->catfile($tmp, "$index.$sample.bedpe");
	
  open (my $ifh, $file) or die "Could not open file '$file' $!";
  open(my $ofh1, '>', $output1) or die "Could not open file '$output1' $!";
  
  while (<$ifh>) {
    chomp;
    my $line = $_;
		  
    my @fields;
    my $name;
    my $chr1;
    my $pos1_start;
    my $pos1_end;
    my $strand1;
    my $chr2;
    my $pos2_start;
    my $pos2_end;
    my $strand2; 
		  
    if($filetype eq 'star'){
      next if($line =~ m/$STAR_HEADER_PATTERN/);
		  	
      @fields = split $STAR_SPLIT_CHAR, $line;
      $name = $fields[$STAR_BREAKREF - 1];
      $chr1 = $fields[$STAR_CHR1 - 1];
      $pos1_start = $fields[$STAR_POS1 - 1]-1;
      $pos1_end = $fields[$STAR_POS1 - 1];
      $strand1 = $fields[$STAR_STRAND1 - 1];
      $chr2 = $fields[$STAR_CHR2 - 1];
      $pos2_start = $fields[$STAR_POS2 - 1]-1;
      $pos2_end = $fields[$STAR_POS2- 1];
      $strand2 = $fields[$STAR_STRAND2 - 1];
		    		    
      print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$chr2."\t".$pos2_start."\t".$pos2_end."\t".$filetype."\t".$name."\t".$strand1."\t".$strand2."\n";
    }
    elsif($filetype eq 'tophat'){
      next if($line =~ m/$TOPHAT_HEADER_PATTERN/);
		    
      @fields = split $TOPHAT_SPLIT_CHAR, $line;
      $name = $fields[$TOPHAT_BREAKREF - 1];
      $chr1 = $fields[$TOPHAT_CHR1 - 1];
      $pos1_start = $fields[$TOPHAT_POS1 - 1]-1;
      $pos1_end = $fields[$TOPHAT_POS1 - 1];
      $strand1 = $fields[$TOPHAT_STRAND1 - 1];
      $chr2 = $fields[$TOPHAT_CHR2 - 1];
      $pos2_start = $fields[$TOPHAT_POS2 - 1]-1;
      $pos2_end = $fields[$TOPHAT_POS2- 1];
      $strand2 = $fields[$TOPHAT_STRAND2 - 1];
		    		    
      print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$chr2."\t".$pos2_start."\t".$pos2_end."\t".$filetype."\t".$name."\t".$strand1."\t".$strand2."\n";
    }
    
    # It must be defuse format
    else{
      next if($line =~ m/$DEFUSE_HEADER_PATTERN/);
		  	
      @fields = split $DEFUSE_SPLIT_CHAR, $line;
      $name = $fields[$DEFUSE_BREAKREF - 1]."_".$fields[$DEFUSE_CLUSTER_ID - 1];
      $chr1 = $fields[$DEFUSE_CHR1 - 1];
      $pos1_start = $fields[$DEFUSE_POS1 - 1]-1;
      $pos1_end = $fields[$DEFUSE_POS1 - 1];
      $strand1 = $fields[$DEFUSE_STRAND1 - 1];
      $chr2 = $fields[$DEFUSE_CHR2 - 1];
      $pos2_start = $fields[$DEFUSE_POS2 - 1]-1;
      $pos2_end = $fields[$DEFUSE_POS2- 1];
      $strand2 = $fields[$DEFUSE_STRAND2 - 1];
		    		    
      print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$chr2."\t".$pos2_start."\t".$pos2_end."\t".$filetype."\t".$name."\t".$strand1."\t".$strand2."\n";
    }
  }
  close ($ifh);
  close ($ofh1);
	
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
	
  return 1;  
  
}

sub deduplicate_fusions {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  my $sample = $options->{'sample'};
  
  my $input_file = File::Spec->catfile($tmp, "$sample.final");
  my $temp_file = File::Spec->catfile($tmp, "$sample.final.2");
  
  # Open the file, read through and add a sort key based on the source algorithm. STD should be 1 and everything else 2.  
  open(my $ifh1, $input_file) or die "Could not open file $input_file $!";
  open(my $ofh1, '>', $temp_file) or die "Could not open file $temp_file $!";
  while(<$ifh1>){
    chomp;
    my $line = $_;
    if($line =~ m/STD$/){
      print $ofh1 "1\t".$line."\n";
    }
    else{
      print $ofh1 "2\t".$line."\n";
    }
  }
  close($ofh1);
  close($ifh1);
  
  my $sorted_file = File::Spec->catfile($tmp, "$sample.final.2.sorted");
  my $sort_command = sprintf $TRI_SORT, 1, 1,8,8,15,15, $temp_file, $sorted_file;
  system($sort_command);
  
  my %seen;
  my $output_file = File::Spec->catfile($tmp, "$sample.final.deduped");
  
  # Read through the sorted file and check there are no duplicate tophat breakpoints (sometimes star and tophat breakpoints differ). 
  open(my $ifh2, $sorted_file) or die "Could not open file $sorted_file $!";
  open(my $ofh2, '>', $output_file) or die "Could not open file $output_file $!";
  while(<$ifh2>){
    chomp;
    my $line = $_;
    my @fields = split "\t", $line;
    my $length = scalar @fields;
    my $sort_key = $fields[0];
    my $source = $fields[$length-1];
    if(defined $sort_key){
      if($sort_key == 1 ){
        $seen{$fields[2]} = $fields[3];
        shift @fields;
        print $ofh2 join("\t", @fields)."\n";
      }
      elsif($source eq 'ST'){
        if(!exists $seen{$fields[2]}){
          $seen{$fields[2]} = $fields[3];
          shift @fields;
          print $ofh2 join("\t", @fields)."\n";
        }
      }     
      elsif($source eq 'D'){
        $fields[1] =~ m/^(.*)_[0-9]+$/;
        my $junction = $1;
        if(!exists $seen{$junction}){
          $seen{$junction} = $fields[1];
          shift @fields;
          print $ofh2 join("\t", @fields)."\n";
        }
      }
      else{
        if(!exists $seen{$fields[1]}){
          $seen{$fields[1]} = $fields[2];
          shift @fields;
          print $ofh2 join("\t", @fields)."\n";
        }
      }
    }
  }
  close($ofh2);
  close($ifh2);
  
  return 1;
}

sub filter_gtf {
  my ($options, $feature) = @_;
	
	my $tmp = $options->{'tmp'};
	
	my $gtf = $options->{'gtf'};
	
  my $filtered_gtf = File::Spec->catfile($tmp, "filtered_$feature.gtf");
  
  unless (-e $filtered_gtf){
  open (my $ifh, $gtf) or die "Could not open file '$gtf' $!";
  open(my $ofh, '>', $filtered_gtf) or die "Could not open file '$filtered_gtf' $!";

  while (<$ifh>) {
    chomp;
    my $line = $_;
    $line =~ s/"//g;
    $line =~ s/;$//;
    next if($line =~ m/^#/ );
		
    my %annotation;
    my @fields = split '\t', $line;
    next if($fields[2] ne $feature);
    my $annot_column = scalar @fields;
    my @annot_fields = split /; /, $fields[$annot_column-1];
    foreach my $item(@annot_fields) {
      my ($type,$value)= split / /, $item;
      $annotation{$type} = $value;
    }
    #print $ofh $line."\n" if(exists $ALLOWED_BIOTYPES{$annotation{'gene_biotype'}}); # UNCOMMENT THIS LINE TO FILTER ON BIOTYPE AND COMMENT OUT LINE BELOE   
    print $ofh $line."\n"
  }
  close($ifh);
  close($ofh);
  }
  return $filtered_gtf;
}

sub find_closest_boundary {
  my ($break, $start, $end) = @_;
  
  my $distance = abs($start - $break);
  my $distance2 = abs($end - $break);
  
  $distance = $distance2 if($distance2 < $distance);
  
  return $distance;
}

sub generate_output {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $sample = $options->{'sample'};

  my $star_file = $options->{'fusion_files'}->{'1'}->{'name'};
  my $tophat_file = $options->{'fusion_files'}->{'2'}->{'name'};
  my $defuse_file = $options->{'fusion_files'}->{'3'}->{'name'};

  my $star_data = parse_star_file($star_file);
  my $tophat_data = parse_tophat_file($tophat_file);
  my $defuse_data = parse_defuse_file($defuse_file);
  
  my $annot_file = File::Spec->catfile($tmp, "$sample.final.deduped");
  my $output_file = File::Spec->catfile($tmp, "$sample.detected.fusions.txt");
  open(my $ofh1, '>', $output_file) or die "Could not open file $output_file $!";
  
  if(-s $annot_file){
  
  open(my $ifh1, $annot_file) or die "Could not open file $annot_file $!";
  print $ofh1 $OUTPUT_HEADER;
  while(<$ifh1>){
    chomp;
    my $line = $_;
    my @fields = split "\t", $line;
    
  	my $length = scalar @fields;
  	my $source = $fields[$length-1];
  	my $confidence = $CONFIDENCE_SCORES{$source};
  	
  	my $star_pos = index($source,'S') if($source =~ m/S/);
  	my $tophat_pos = index($source,'T') if($source =~ m/T/);
  	my $defuse_pos = index($source,'D') if($source =~ m/D/);
  	
  	my $star_breakpoint = 'NA';
  	my $tophat_breakpoint = 'NA';
  	my $defuse_breakpoint = 'NA';
  	my $defuse_junction = 'NA';
  	my $defuse_clusterid = 'NA';	
  	my $chr1 = 'NA';
  	my $pos1 = 'NA';
  	my $strand1 = 'NA';
  	my $chr2 = 'NA';
  	my $pos2 = 'NA';
  	my $strand2 = 'NA';
   	my $star_junction_reads = 'NA';
  	my $star_spanning_frags = 'NA';
  	my $tophat_junction_reads = 'NA';
  	my $tophat_spanning_frags = 'NA';
  	my $tophat_splitr_seq = 'NA';
  	my $defuse_splitr_count = 'NA';
  	my $defuse_span_count = 'NA';
  	my $defuse_splitr_seq = 'NA';
  	  
  	if(defined $tophat_pos){
  	  $tophat_breakpoint = $fields[$tophat_pos];
  	  $tophat_junction_reads = $tophat_data->{$tophat_breakpoint}{'num_spanning_reads'};
  	  $tophat_spanning_frags = $tophat_data->{$tophat_breakpoint}{'num_spanning_mate_pairs'};
  	  $chr1 = $tophat_data->{$tophat_breakpoint}{'chr1'};
  	  $pos1 = $tophat_data->{$tophat_breakpoint}{'pos1'};
  	  $strand1 = $tophat_data->{$tophat_breakpoint}{'strand1'};
  	  $chr2 = $tophat_data->{$tophat_breakpoint}{'chr2'};
  	  $pos2 = $tophat_data->{$tophat_breakpoint}{'pos2'};
  	  $strand2 = $tophat_data->{$tophat_breakpoint}{'strand2'};
  	} 	  
  	if(defined $star_pos){
  	  $star_breakpoint = $fields[$star_pos];  
  	  $star_junction_reads = $star_data->{$star_breakpoint}{'junction_reads'};
  	  $star_spanning_frags = $star_data->{$star_breakpoint}{'spanning_frags'};
  	  $chr1 = $star_data->{$star_breakpoint}{'chr1'};
  	  $pos1 = $star_data->{$star_breakpoint}{'pos1'};
  	  $strand1 = $star_data->{$star_breakpoint}{'strand1'};
  	  $chr2 = $star_data->{$star_breakpoint}{'chr2'};
  	  $pos2 = $star_data->{$star_breakpoint}{'pos2'};
  	  $strand2 = $star_data->{$star_breakpoint}{'strand2'};
  	}
  	if(defined $defuse_pos){
  	  $defuse_breakpoint = $fields[$defuse_pos];
  	  $defuse_splitr_count = $defuse_data->{$defuse_breakpoint}{'split_reads'};
  	  $defuse_span_count = $defuse_data->{$defuse_breakpoint}{'span_reads'};
  	  # If we already have the chr-pos-strand info we don't want to pick it up from deFuse in case the orientation is reported differently
  	  if(!defined $star_pos && !defined $tophat_pos){
 	      $chr1 = $defuse_data->{$defuse_breakpoint}{'chr1'};
  	    $pos1 = $defuse_data->{$defuse_breakpoint}{'pos1'};
  	    $strand1 = $defuse_data->{$defuse_breakpoint}{'strand1'};
  	    $chr2 = $defuse_data->{$defuse_breakpoint}{'chr2'};
  	    $pos2 = $defuse_data->{$defuse_breakpoint}{'pos2'};
  	    $strand2 = $defuse_data->{$defuse_breakpoint}{'strand2'};
  	  }
  	  $defuse_splitr_seq = $defuse_data->{$defuse_breakpoint}{'sequence'};
  	  my @defuse_temp = split "_", $defuse_breakpoint;
  	  $defuse_junction = $defuse_temp[0];
  	  $defuse_clusterid = $defuse_temp[1];
  	}
    if($length > 11){
  	  my $gene1_name = $fields[3];
  	  my $gene1_id = $fields[4];
  	  my $gene2_name = $fields[10];
  	  my $gene2_id = $fields[11];
      my $fusion_name = $gene1_name."--".$gene2_name;
      my $transcript1_id = $fields[5];
      my $transcript1_src = $fields[6];
      my $exon1_number = $fields[7];
      my $exon1_start = $fields[8];
      my $exon1_end = $fields[9];
      my $transcript2_id = $fields[12];
      my $transcript2_src = $fields[13];
      my $exon2_number = $fields[14];
      my $exon2_start = $fields[15];
      my $exon2_end = $fields[16];
      
      $source = reverse $source;
      
      print $ofh1 "$sample\t$fusion_name\t$source\t$confidence\t$star_breakpoint\t$tophat_breakpoint\t$defuse_junction\t$defuse_clusterid\t$star_junction_reads\t$star_spanning_frags\t$tophat_junction_reads\t$tophat_spanning_frags\t$defuse_splitr_count\t$defuse_span_count\t$gene1_name\t$gene1_id\t$chr1\t$pos1\t$strand1\t$gene2_name\t$gene2_id\t$chr2\t$pos2\t$strand2\t$transcript1_id\t$transcript1_src\t$exon1_number\t$exon1_start\t$exon1_end\t$transcript2_id\t$transcript2_src\t$exon2_number\t$exon2_start\t$exon2_end\t$defuse_splitr_seq\t$tophat_splitr_seq\n";
    }
    else{
      $source = reverse $source;
      print $ofh1 "$sample\tFUSION COULD NOT BE ANNOTATED\t$source\t$confidence\t$star_breakpoint\t$tophat_breakpoint\t$defuse_junction\t$defuse_clusterid\t$star_junction_reads\t$star_spanning_frags\t$tophat_junction_reads\t$tophat_spanning_frags\t$defuse_splitr_count\t$defuse_span_count\t\t\t$chr1\t$pos1\t$strand1\t\t\t$chr2\t$pos2\t$strand2\t\t\t\t\t\t\t\t\t\t\t$defuse_splitr_seq\t$tophat_splitr_seq\n";
    }
  }
  close($ifh1);
  }
  
  close($ofh1);
  
  return 1;
}

sub parse_annotation {
  my $line = shift;
  $line =~ s/"//g;
  
  my %annotation;
  my @fields = split "\t", $line;
  
  $annotation{'breakpoint'} = $fields[3];
  $annotation{'alt_breakpoint'} = $fields[4];
  $annotation{'alt_breakpoint2'} = $fields[8];
  $annotation{'chr'} = $fields[0];
  $annotation{'pos_start'} = $fields[1];
  $annotation{'pos_end'} = $fields[2];
  $annotation{'strand'} = $fields[5];
  $annotation{'feature'} = $fields[12];
  $annotation{'feature_start'} = $fields[13];
  $annotation{'feature_end'} = $fields[14];
  $annotation{'genename'} = $fields[6];
  $annotation{'source'} = $fields[9];

	my $annot_column = scalar @fields;
  my @annot_fields = split /; /, $fields[$annot_column-1];
  foreach my $item(@annot_fields) {
    my ($type,$value)= split / /, $item;
    $annotation{$type} = $value;
  }

  return \%annotation;
}

sub parse_bed_file {
  my ($line, $breaknum) = @_;
  
  my @fields = split "\t", $line; 
  my $fusion;
  
  if($breaknum == 1){
  
    $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
    -breakpoint	=> $fields[3],
    -alt_breakpoint	=> $fields[4],
    -alt_breakpoint2	=> $fields[8],
    -chr1	=> $fields[0],
    -pos1_start	=> $fields[1],
    -pos1_end	=> $fields[2],
    -strand1	=> $fields[5],
    -gene1	=> $fields[6],
    -gene1_id	=> $fields[7],
    -feature1 => $fields[9]);
  }
  else{
    $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
    -breakpoint	=> $fields[3],
    -alt_breakpoint	=> $fields[4],
    -alt_breakpoint2	=> $fields[8],
    -chr2	=> $fields[0],
    -pos2_start	=> $fields[1],
    -pos2_end	=> $fields[2],
    -strand2	=> $fields[5],
    -gene2	=> $fields[6],
    -gene2_id	=> $fields[7],
    -feature1 => $fields[9]);    
  }
  
  return $fusion;
}

sub parse_break_data {
  my $line = shift;
  
  my %break;
  my @fields = split "\t", $line;
  
  my @breakpoint = split "_", $fields[0];
  $break{'breakpoint'} = $breakpoint[0];
  $break{'alt_breakpoint'} = $fields[1];
  $break{'gene_name'} = $fields[2];
  $break{'gene_id'} = $fields[3];
  $break{'strand'} = $fields[4];
  $break{'feature'} = $fields[5];
  $break{'exon_id'} = $fields[6];
  $break{'exon_number'} = $fields[7];
  $break{'exon_start'} = $fields[8];
  $break{'exon_end'} = $fields[9];
  $break{'transcript_id'} = $fields[10];
  $break{'transcript_src'} = $fields[11];
  $break{'gene_biotype'} = $fields[12];
  
  return \%break;
}

sub parse_defuse_file {
  
  my $defuse_file = shift;
  
  my %defuse_data;
  open (my $ifh1, $defuse_file) or die "Could not open file '$defuse_file' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    next if($line =~ m/$DEFUSE_HEADER_PATTERN/);
    my @fields = split $DEFUSE_SPLIT_CHAR, $line;
    my $break_ref = $fields[$DEFUSE_BREAKREF-1];
    my $cluster_id = $fields[$DEFUSE_CLUSTER_ID-1];
    my $breakpoint = $break_ref."_".$cluster_id;
    $defuse_data{$breakpoint}{'breakpoint'} = $break_ref;
    $defuse_data{$breakpoint}{'cluster_id'} = $cluster_id;
    $defuse_data{$breakpoint}{'chr1'} = $fields[$DEFUSE_CHR1-1];
    $defuse_data{$breakpoint}{'pos1'} = $fields[$DEFUSE_POS1-1];
    $defuse_data{$breakpoint}{'strand1'} = $fields[$DEFUSE_STRAND1-1];
    $defuse_data{$breakpoint}{'gene1_name'} = $fields[$DEFUSE_GENENAME1-1];
    $defuse_data{$breakpoint}{'gene1_id'} = $fields[$DEFUSE_GENEID1-1];
    $defuse_data{$breakpoint}{'chr2'} = $fields[$DEFUSE_CHR2-1];
    $defuse_data{$breakpoint}{'pos2'} = $fields[$DEFUSE_POS2-1];
    $defuse_data{$breakpoint}{'strand2'} = $fields[$DEFUSE_STRAND2-1];
    $defuse_data{$breakpoint}{'gene2_name'} = $fields[$DEFUSE_GENENAME2-1];
    $defuse_data{$breakpoint}{'gene2_id'} = $fields[$DEFUSE_GENEID2-1];
    $defuse_data{$breakpoint}{'sequence'} = $fields[$DEFUSE_SEQUENCE-1];
    $defuse_data{$breakpoint}{'split_reads'} = $fields[$DEFUSE_SPLIT_READS-1];
    $defuse_data{$breakpoint}{'span_reads'} = $fields[$DEFUSE_SPAN_READS-1];
  }
  close ($ifh1);

  return \%defuse_data;
}

sub parse_exon_data {
  my ($fusion, $breaknum, $exons) = @_;

	my $exon_number;
	my $exon_start;
	my $exon_end;
	my $transcript_id;
	my $gene_biotype;
	my $found_exon_boundary = 0;
	my $found_exon = 0;
	my $curr_distance = 10000000;

  my $break_pos = $fusion->{'pos'.$breaknum.'_end'};

  for my $e (keys $exons){
    my $exon = $exons->{$e};
    if($break_pos == $exon->{'feature_start'} || $break_pos == $exon->{'feature_end'}){
      $found_exon_boundary = 1;
      $transcript_id = $exon->{'transcript_id'};
		  $gene_biotype = $exon->{'gene_biotype'};
		  $exon_start = $exon->{'feature_start'};
		  $exon_end = $exon->{'feature_end'};
		  $exon_number = $exon->{'exon_number'};
      last;
    }
    elsif($break_pos > $exon->{'feature_start'} && $break_pos < $exon->{'feature_end'}){
      my $distance = find_closest_boundary($break_pos, $exon->{'feature_start'}, $exon->{'feature_end'});
      if ($distance < $curr_distance){
		    $found_exon = 1;
        $transcript_id = $exon->{'transcript_id'};
		    $gene_biotype = $exon->{'gene_biotype'};
		    $exon_start = $exon->{'feature_start'};
		    $exon_end = $exon->{'feature_end'};
		    $exon_number = $exon->{'exon_number'};
		  }
    }
    last if($found_exon_boundary);
  }
	if ($found_exon_boundary || $found_exon){
	  if($breaknum == 1){
	    $fusion->transcript1_id($transcript_id);
		  $fusion->gene1_biotype($gene_biotype);
		  $fusion->exon1_num($exon_number);
		  $fusion->feature1_start($exon_start);
		  $fusion->feature1_end($exon_end);
	  }
	  else{
	    $fusion->transcript2_id($transcript_id);
		  $fusion->gene2_biotype($gene_biotype);
		  $fusion->exon2_num($exon_number);
		  $fusion->feature2_start($exon_start);
		  $fusion->feature2_end($exon_end);
    }
  }
  return $fusion;
}

sub parse_gene_info {
  my $line = shift;
  $line =~ s/"//g;
  
  my %gene_annotation;
  my @fields = split "\t", $line;
  
  $gene_annotation{'start'} = $fields[3];
  $gene_annotation{'end'} = $fields[4];
  
  my $annot_column = scalar @fields;
  my @annot_fields = split /; /, $fields[$annot_column-1];
  foreach my $item(@annot_fields) {
    my ($type,$value)= split / /, $item;
    $gene_annotation{$type} = $value;
  }
  
  return \%gene_annotation;
}

sub parse_intersection {
  my ($line) = @_;
  
  my @fields = split " ", $line;
  
  my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation();
  
  $fusion->breakpoint($fields[0]);
  $fusion->chr1($fields[1]);
  $fusion->pos1_start($fields[2]);
  $fusion->pos1_end($fields[3]);
  $fusion->strand1($fields[8]);
  $fusion->chr2($fields[4]);
  $fusion->pos2_start($fields[5]);
  $fusion->pos2_end($fields[6]);
  $fusion->strand2($fields[9]);
  $fusion->alt_breakpoint($fields[17]);
  $fusion->alt_breakpoint2($fields[36]);
  
  return $fusion;
}

sub parse_overlap {
  my ($line) = @_;
  
  my @fields = split "\t", $line;
  
  my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
    -breakpoint	=> $fields[7],
    -chr1	=> $fields[0],
    -pos1_start	=> $fields[1],
    -pos1_end	=> $fields[2],
    -strand1	=> $fields[8],
    -chr2	=> $fields[3],
    -pos2_start	=> $fields[4],
    -pos2_end	=> $fields[5],
    -strand2	=> $fields[9],
    -alt_breakpoint => $fields[17],
    -alt_breakpoint2 => 'NA');
  
  return $fusion;
}

sub parse_star_file {
  
  my $star_file = shift;
  
  my %star_data;
  open (my $ifh1, $star_file) or die "Could not open file '$star_file' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    next if($line =~ m/$STAR_HEADER_PATTERN/);
    my @fields = split $STAR_SPLIT_CHAR, $line;
    my $breakpoint = $fields[0];
    $star_data{$breakpoint}{'fusion_name'} = $fields[$STAR_FUSION_NAME-1];
    $star_data{$breakpoint}{'junction_reads'} = $fields[$STAR_JUNCTION_READS-1];
    $star_data{$breakpoint}{'spanning_frags'} = $fields[$STAR_SPANNING_FRAGS-1];
    $star_data{$breakpoint}{'gene1_name'} = $fields[$STAR_GENENAME1-1];
    $star_data{$breakpoint}{'gene1_id'} = $fields[$STAR_GENEID1-1];
    $star_data{$breakpoint}{'chr1'} = $fields[$STAR_CHR1-1];
    $star_data{$breakpoint}{'pos1'} = $fields[$STAR_POS1-1];
    $star_data{$breakpoint}{'strand1'} = $fields[$STAR_STRAND1-1];
    $star_data{$breakpoint}{'dis_exon_1'} = $fields[$STAR_DIS_EXON1-1];
    $star_data{$breakpoint}{'gene2_name'} = $fields[$STAR_GENENAME2-1];
    $star_data{$breakpoint}{'gene2_id'} = $fields[$STAR_GENEID2-1];
    $star_data{$breakpoint}{'chr2'} = $fields[$STAR_CHR2-1];
    $star_data{$breakpoint}{'pos2'} = $fields[$STAR_POS2-1];
    $star_data{$breakpoint}{'strand2'} = $fields[$STAR_STRAND2-1];
    $star_data{$breakpoint}{'dis_exon_2'} = $fields[$STAR_DIS_EXON2-1];
  }
  close ($ifh1);

  return \%star_data;
}

sub parse_tophat_file {
  
  my $tophat_file = shift;
  
  my %tophat_data;
  open (my $ifh1, $tophat_file) or die "Could not open file '$tophat_file' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    next if($line =~ m/$TOPHAT_HEADER_PATTERN/);
    my @fields = split $TOPHAT_SPLIT_CHAR, $line;
    my $breakpoint = $fields[$TOPHAT_BREAKREF-1];
    $tophat_data{$breakpoint}{'gene1_name'} = $fields[$TOPHAT_GENE1-1];
    $tophat_data{$breakpoint}{'chr1'} = $fields[$TOPHAT_CHR1-1];
    $tophat_data{$breakpoint}{'pos1'} = $fields[$TOPHAT_POS1-1];
    $tophat_data{$breakpoint}{'gene2_name'} = $fields[$TOPHAT_GENE2-1];
    $tophat_data{$breakpoint}{'chr2'} = $fields[$TOPHAT_CHR2-1];
    $tophat_data{$breakpoint}{'pos2'} = $fields[$TOPHAT_POS2-1];
    $tophat_data{$breakpoint}{'num_spanning_reads'} = $fields[$TOPHAT_SPAN_READS-1];
    $tophat_data{$breakpoint}{'num_spanning_mate_pairs'} = $fields[$TOPHAT_SPAN_MATE_PAIRS-1];
    $tophat_data{$breakpoint}{'num_spanning_mates_2'} = $fields[$TOPHAT_SPAN_MATE_PAIRS2-1];
    $tophat_data{$breakpoint}{'score'} = $fields[$TOPHAT_SCORE];
    $tophat_data{$breakpoint}{'strand1'} = $fields[$TOPHAT_STRAND1-1];
    $tophat_data{$breakpoint}{'strand2'} = $fields[$TOPHAT_STRAND2-1];
  }
  close ($ifh1);

  return \%tophat_data;
}

sub parse_transcript_data {
  my ($fusion, $breaknum, $transcripts) = @_;
  my @filteredTrans;
  foreach my $t(@{$transcripts}){
		push(@filteredTrans, $t) if ($fusion->{'pos'.$breaknum.'_end'} >= $t->getGenomicMinPos && $fusion->{'pos'.$breaknum.'_end'} <= $t->getGenomicMaxPos);
	}
	my @sortedTrans = sort{&annotation_sort} @filteredTrans;
	
	if(defined $sortedTrans[0]){
    my $vagrent_genename;
	  my $exon_number;
		my $exon_start;
		my $exon_end;
		my $transcript_id;
		my $gene_biotype;
		my $found_exon_boundary = 0;
		my $found_exon = 0;
		my $curr_distance = 10000000;
		my $num_transcripts = scalar @sortedTrans;
		for (my $x=0;$x<$num_transcripts; $x++){
		  my @exons = $sortedTrans[$x]->getExons;
		  my $num_exons = scalar @exons;
		  for (my $y=0;$y<$num_exons; $y++){
		    my $e = $exons[$y];
		    if($fusion->{'pos'.$breaknum.'_end'} == $e->getMinPos || $fusion->{'pos'.$breaknum.'_end'} == $e->getMaxPos){
		      $found_exon_boundary = 1;
		      $transcript_id = $sortedTrans[$x]->getAccession;
		      $gene_biotype = $sortedTrans[$x]->{'_genetype'};
		      $exon_start = $e->getMinPos;
		      $exon_end = $e->getMaxPos;
		      $exon_number = $y+1;
					$vagrent_genename = $sortedTrans[$x]->{'_genename'};
		      last;
		    }elsif($fusion->{'pos'.$breaknum.'_end'} > $e->getMinPos && $fusion->{'pos'.$breaknum.'_end'} < $e->getMaxPos){
		       my $distance = find_closest_boundary($fusion->{'pos'.$breaknum.'_end'}, $e->getMinPos, $e->getMaxPos);
		       if ($distance < $curr_distance){
		         $found_exon = 1;
		         $transcript_id = $sortedTrans[$x]->getAccession;
		         $gene_biotype = $sortedTrans[$x]->{'_genetype'};
		         $exon_start = $e->getMinPos;
		         $exon_end = $e->getMaxPos;
		         $exon_number = $y+1;
		  			 $vagrent_genename = $sortedTrans[$x]->{'_genename'};
		       }
		    }
		  }
		  last if($found_exon_boundary);
		}

  if(defined $vagrent_genename){
    if($vagrent_genename eq $fusion->{'gene'.$breaknum} && ($found_exon_boundary || $found_exon)){
		    if($breaknum == 1){
		      $fusion->transcript1_id($transcript_id);
		      $fusion->gene1_biotype($gene_biotype);
		      $fusion->exon1_num($exon_number);
		      $fusion->feature1_start($exon_start);
		      $fusion->feature1_end($exon_end);
		    }
		    else{
			    $fusion->transcript2_id($transcript_id);
		      $fusion->gene2_biotype($gene_biotype);
		      $fusion->exon2_num($exon_number);
		      $fusion->feature2_start($exon_start);
		      $fusion->feature2_end($exon_end);
        }
	    }
    }	
  }
  return $fusion;
}

sub parse_vagrent_query_file {
  my ($line) = @_;
  
  my @fields = split "\t", $line;
  
  my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
    -breakpoint	=> $fields[0],
    -alt_breakpoint	=> $fields[1],
    -alt_breakpoint2	=> $fields[2],
    -chr1	=> $fields[3],
    -pos1_start	=> $fields[4],
    -pos1_end	=> $fields[5],
    -strand1	=> $fields[6],
    -chr2	=> $fields[7],
    -pos2_start	=> $fields[8],
    -pos2_end	=> $fields[9],
    -strand2	=> $fields[10],
    -feature1 => $fields[11]);
  
  return $fusion;
}

sub process_annotation_file {
  my ($input) = @_;

  my %exon_annotation;
  
  open (my $ifh1, $input) or die "Could not open file '$input' $!";
  while (<$ifh1>){
    chomp;
    my $line = $_;
    my $annotation = parse_annotation($line);
    if($annotation->{'gene_name'} eq $annotation->{'genename'}){
      my $break = $annotation->{'breakpoint'};
      my $exon = $annotation->{'exon_id'};
      $exon_annotation{$break}{$exon} = $annotation;
    }
  }
  close ($ifh1);

  return \%exon_annotation;
}

sub process_overlap_files {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $sample = $options->{'sample'};
  
  # Establish the source of 1 and 2 respectively
  my $source_comb;
  my $source1 = $options->{'fusion_files'}->{'1'}->{'format'};
  my $source2 = $options->{'fusion_files'}->{'2'}->{'format'};
  
  my $output_file = File::Spec->catfile($tmp, "$sample.vagrent.query.list");
  my %all_fusions;
  my $cols;
  
  if($options->{'num'} == 3){
    my $source3 = $options->{'fusion_files'}->{'3'}->{'format'};
    $source_comb = uc(substr($source1,0,1).substr($source2,0,1).substr($source3,0,1));
    my $overlap_file1_2_3 = File::Spec->catfile($tmp, "1_2_3.$sample.bedpe_overlap");
    open (my $ifh1, $overlap_file1_2_3) or die "Could not open file '$overlap_file1_2_3' $!";
    while (<$ifh1>) {
      chomp;
      my $line = $_;
      my $fusion = parse_intersection($line);
      my $breakpoint = $fusion->breakpoint();
      
      if(!exists $all_fusions{$breakpoint}){
        $all_fusions{$breakpoint} = $fusion;
        $all_fusions{$breakpoint}{'source'} = $source_comb;
      }
    }
    close($ifh1);
    
    $source_comb = uc(substr($source1,0,1).substr($source3,0,1));
    my $overlap_file1_3 = File::Spec->catfile($tmp, "1_3.$sample.bedpe_overlap");
    open (my $ifh2, $overlap_file1_3) or die "Could not open file '$overlap_file1_3' $!";
    while (<$ifh2>) {
      chomp;
      my $line = $_;
      my $fusion = parse_overlap($line);
      my $breakpoint = $fusion->breakpoint();
      
      if(!exists $all_fusions{$breakpoint}){
        $all_fusions{$breakpoint} = $fusion;
        $all_fusions{$breakpoint}{'source'} = $source_comb;
      }
    }
    close($ifh2); 
    
    $source_comb = uc(substr($source2,0,1).substr($source3,0,1));
    my $overlap_file2_3 = File::Spec->catfile($tmp, "2_3.$sample.bedpe_overlap");
    open (my $ifh3, $overlap_file2_3) or die "Could not open file '$overlap_file2_3' $!";
    while (<$ifh3>) {
      chomp;
      my $line = $_;
      my $fusion = parse_overlap($line);
      my $breakpoint = $fusion->breakpoint();
      
      if(!exists $all_fusions{$breakpoint}){
        $all_fusions{$breakpoint} = $fusion;
        $all_fusions{$breakpoint}{'source'} = $source_comb;
      }
    }
    close($ifh3);    
  }
  
  my $overlap_file1_2 = File::Spec->catfile($tmp, "1_2.$sample.bedpe_overlap");
  $source_comb = uc(substr($source1,0,1).substr($source2,0,1));
  open (my $ifh4, $overlap_file1_2) or die "Could not open file '$overlap_file1_2' $!";
  while (<$ifh4>) {
    chomp;
    my $line = $_;
    my $fusion = parse_overlap($line);
    my $breakpoint = $fusion->breakpoint();
      
    if(!exists $all_fusions{$breakpoint}){
      $all_fusions{$breakpoint} = $fusion;
      $all_fusions{$breakpoint}{'source'} = $source_comb;
    }
  }
  close($ifh4);
  
  open(my $ofh1, '>', $output_file) or die "Could not open file '$output_file' $!";
  for my $brk (keys %all_fusions){
    my $output_line = $all_fusions{$brk}->format_fusion_line($all_fusions{$brk}{'source'});
    print $ofh1 $output_line."\n";
  
  }
  close($ofh1);
  
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  return 1;
}

sub process_singletons {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  my $sample = $options->{'sample'};
  
  my %star_ids;
  my %tophat_ids;
  my %defuse_ids;
  
  my $overlaps_file = File::Spec->catfile($tmp, "$sample.vagrent.query.list");
  open (my $ifh1, $overlaps_file) or die "Could not open file '$overlaps_file' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    my @fields = split "\t", $line;
    my $length = scalar @fields;
    my $source = $fields[$length-1];
    
    if($source eq 'STD'){
      $star_ids{$fields[0]} = 'STD';
      $tophat_ids{$fields[1]} = 'STD';
      my $defuse_brk_id = $fields[2];    
      my @defuse_id_split = split "_", $defuse_brk_id;
      my $defuse_junction = $defuse_id_split[0];
      $defuse_ids{$defuse_junction} = 'STD';
    }
    elsif($source eq 'SD'){
      $star_ids{$fields[0]} = 'SD';
      my $defuse_brk_id = $fields[1];    
      my @defuse_id_split = split "_", $defuse_brk_id;
      my $defuse_junction = $defuse_id_split[0];
      $defuse_ids{$defuse_junction} = 'SD';
    }
    elsif($source eq 'ST'){
      $star_ids{$fields[0]} = 'ST';
      $tophat_ids{$fields[1]} = 'ST';
    }
    elsif($source eq 'TD'){
      $tophat_ids{$fields[0]} = 'TD';
      $defuse_ids{$fields[1]} = 'TD';
      my $defuse_brk_id = $fields[1];    
      my @defuse_id_split = split "_", $defuse_brk_id;
      my $defuse_junction = $defuse_id_split[0];
      $defuse_ids{$defuse_junction} = 'TD';
    }
    else{
      die "Fusion does not have recognised source algorithms\n"; 
    }
  }
  close($ifh1);

  my $star_file = $options->{'fusion_files'}->{'1'}->{'name'};
  my $tophat_file = $options->{'fusion_files'}->{'2'}->{'name'};
  my $defuse_file = $options->{'fusion_files'}->{'3'}->{'name'};
  my $star_data = parse_star_file($star_file);
  my $tophat_data = parse_tophat_file($tophat_file);
  my $defuse_data = parse_defuse_file($defuse_file);
  my %all_singletons;
  
  for my $star_brk (keys $star_data){
    if(!exists $star_ids{$star_brk}){
    
      $star_ids{$star_brk} = 'S';
    
      my $pos_start1 = $star_data->{$star_brk}{'pos1'} -1;
      my $pos_start2 = $star_data->{$star_brk}{'pos2'} -1;
    
      my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
      -breakpoint	=> $star_brk,
      -chr1	=> $star_data->{$star_brk}{'chr1'},
      -pos1_start	=> $star_data->{$star_brk}{'pos1'}-1,
      -pos1_end	=> $star_data->{$star_brk}{'pos1'},
      -strand1	=> $star_data->{$star_brk}{'strand1'},
      -chr2	=> $star_data->{$star_brk}{'chr2'},
      -pos2_start	=> $star_data->{$star_brk}{'pos2'}-1,
      -pos2_end	=> $star_data->{$star_brk}{'pos2'},
      -strand2	=> $star_data->{$star_brk}{'strand2'},
      -alt_breakpoint => 'NA',
      -alt_breakpoint2 => 'NA');
      
      if(!exists $all_singletons{$star_brk}){
        $all_singletons{$star_brk} = $fusion;
        $all_singletons{$star_brk}{'source'} = 'S';
      }
    }
  }

  for my $tophat_brk (keys $tophat_data){
    if(!exists $tophat_ids{$tophat_brk}){
    
      $tophat_ids{$tophat_brk} = 'T';
    
      my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
      -breakpoint	=> $tophat_brk,
      -chr1	=> $tophat_data->{$tophat_brk}{'chr1'},
      -pos1_start	=> $tophat_data->{$tophat_brk}{'pos1'}-1,
      -pos1_end	=> $tophat_data->{$tophat_brk}{'pos1'},
      -strand1	=> $tophat_data->{$tophat_brk}{'strand1'},
      -chr2	=> $tophat_data->{$tophat_brk}{'chr2'},
      -pos2_start	=> $tophat_data->{$tophat_brk}{'pos2'}-1,
      -pos2_end	=> $tophat_data->{$tophat_brk}{'pos2'},
      -strand2	=> $tophat_data->{$tophat_brk}{'strand2'},
      -alt_breakpoint => 'NA',
      -alt_breakpoint2 => 'NA');
      
      if(!exists $all_singletons{$tophat_brk}){
        $all_singletons{$tophat_brk} = $fusion;
        $all_singletons{$tophat_brk}{'source'} = 'T';
      }
    }
  }
  
  for my $defuse_brk (keys $defuse_data){
    my @defuse_ids = split "_", $defuse_brk;
    my $breakpoint = $defuse_ids[0];
    if(!exists $defuse_ids{$breakpoint}){
      $defuse_ids{$breakpoint} = 'D';
    
      my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
      -breakpoint	=> $defuse_brk,
      -chr1	=> $defuse_data->{$defuse_brk}{'chr1'},
      -pos1_start	=> $defuse_data->{$defuse_brk}{'pos1'}-1,
      -pos1_end	=> $defuse_data->{$defuse_brk}{'pos1'},
      -strand1	=> $defuse_data->{$defuse_brk}{'strand1'},
      -chr2	=> $defuse_data->{$defuse_brk}{'chr2'},
      -pos2_start	=> $defuse_data->{$defuse_brk}{'pos2'}-1,
      -pos2_end	=> $defuse_data->{$defuse_brk}{'pos2'},
      -strand2	=> $defuse_data->{$defuse_brk}{'strand2'},
      -alt_breakpoint => 'NA',
      -alt_breakpoint2 => 'NA');
      
      if(!exists $all_singletons{$defuse_brk}){
        $all_singletons{$defuse_brk} = $fusion;
        $all_singletons{$defuse_brk}{'source'} = 'D';
      }
    }
  }
  open(my $ofh1, '>>', $overlaps_file) or die "Could not open file '$overlaps_file' $!";
  for my $brk (keys %all_singletons){
    my $output_line = $all_singletons{$brk}->format_fusion_line($all_singletons{$brk}{'source'});
    print $ofh1 $output_line."\n";
  }
  close($ofh1);
  
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  
  return 1;
}

sub query_vagrent {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
  my $sample = $options->{'sample'};

  my $star_file = $options->{'fusion_files'}->{'1'}->{'name'};
  my $tophat_file = $options->{'fusion_files'}->{'2'}->{'name'};
  my $defuse_file = $options->{'fusion_files'}->{'3'}->{'name'};
 
  my $gene_list = parse_star_file($star_file);
  my $gene_gtf = filter_gtf($options, 'gene');

  my %gene_info;
  open (my $ifh1, $gene_gtf) or die "Could not open file '$gene_gtf' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    my $gene_annot = parse_gene_info($line);
    if(!exists $gene_info{$gene_annot->{'gene_name'}}){
      $gene_info{$gene_annot->{'gene_name'}} = $gene_annot->{'gene_id'};
    }
  }
  close ($ifh1);
  
  my %rev_gene_info = reverse %gene_info;
  
  my $tophat_data = parse_tophat_file($tophat_file);
  
  for my $brk (keys $tophat_data){
    if(!exists $gene_list->{$brk}){
      my $gene1_name = $tophat_data->{$brk}{'gene1_name'};
      my $gene2_name = $tophat_data->{$brk}{'gene2_name'};
      $gene_list->{$brk}{'gene1_name'} = $gene1_name;
      if(defined $gene_info{$gene1_name}){
        $gene_list->{$brk}{'gene1_id'} = $gene_info{$gene1_name};
			}
			elsif(defined $rev_gene_info{$gene1_name}){
			  $gene_list->{$brk}{'gene1_id'} = $gene1_name;
			  $gene_list->{$brk}{'gene1_name'} = $rev_gene_info{$gene1_name};
			}
			else{
			  $gene_list->{$brk}{'gene1_id'} = $gene1_name;
			}
      $gene_list->{$brk}{'gene2_name'} = $gene2_name;
      if(defined $gene_info{$gene2_name}){
        $gene_list->{$brk}{'gene2_id'} = $gene_info{$gene2_name};
			}
			elsif(defined $rev_gene_info{$gene2_name}){
			  $gene_list->{$brk}{'gene2_id'} = $gene2_name;
			  $gene_list->{$brk}{'gene2_name'} = $rev_gene_info{$gene2_name};
			}
			else{
			  $gene_list->{$brk}{'gene2_id'} = $gene1_name;
			}
    }
  }
  
  my $defuse_data = parse_defuse_file($defuse_file);
  
  for my $defuse_brk (keys $defuse_data){
    if(!exists $gene_list->{$defuse_brk}){
      my $gene1_name = $defuse_data->{$defuse_brk}{'gene1_name'};
      my $gene2_name = $defuse_data->{$defuse_brk}{'gene2_name'};
      $gene_list->{$defuse_brk}{'gene1_name'} = $gene1_name;
      $gene_list->{$defuse_brk}{'gene1_id'} = $gene_info{$gene1_name};
      $gene_list->{$defuse_brk}{'gene1_id'} = $gene1_name if(!defined $gene_list->{$defuse_brk}{'gene1_id'});
      $gene_list->{$defuse_brk}{'gene2_name'} = $gene2_name;
      $gene_list->{$defuse_brk}{'gene2_id'} = $gene_info{$gene2_name};
      $gene_list->{$defuse_brk}{'gene2_id'} = $gene2_name if(!defined $gene_list->{$defuse_brk}{'gene2_id'});
    }
  }
  
	my $vagrent_version = "VAGrENT_".Sanger::CGP::Vagrent->VERSION;
	
	my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $options->{'cache'});	
	my %breaklist;
	
	my $vagrent_query_file = File::Spec->catfile($tmp, "$sample.vagrent.query.list");

  open (my $ifh2, $vagrent_query_file) or die "Could not open file '$vagrent_query_file' $!";
  while (<$ifh2>) {
    chomp;
    my $line = $_;
    my $fusion = parse_vagrent_query_file($line);
    $fusion->gene1($gene_list->{$fusion->{'breakpoint'}}{'gene1_name'});
    $fusion->gene2($gene_list->{$fusion->{'breakpoint'}}{'gene2_name'});
    $fusion->gene1_id($gene_list->{$fusion->{'breakpoint'}}{'gene1_id'});
    $fusion->gene2_id($gene_list->{$fusion->{'breakpoint'}}{'gene2_id'});
    
    my $genomic_pos1 = Sanger::CGP::Vagrent::Data::GenomicRegion->new('species' => 'human', 'genomeVersion' => 'GRCh38', 'chr' => $fusion->{'chr1'}, 'minpos' => $fusion->{'pos1_start'}, 'maxpos' => $fusion->{'pos1_end'}, 'id' => $fusion->{'breakpoint'});
    my $genomic_pos2 = Sanger::CGP::Vagrent::Data::GenomicRegion->new('species' => 'human', 'genomeVersion' => 'GRCh38', 'chr' => $fusion->{'chr2'}, 'minpos' => $fusion->{'pos2_start'}, 'maxpos' => $fusion->{'pos2_end'}, 'id' => $fusion->{'breakpoint'});
	  
	  unless($genomic_pos1->{'_chr'} eq 'GL000219.1' || $genomic_pos2->{'_chr'} eq 'GL000219.1' || $genomic_pos1->{'_chr'} eq 'KI270726.1' || $genomic_pos2->{'_chr'} eq 'KI270726.1'){
	  
      my @trans1 = $ts->getTranscripts($genomic_pos1);
      my @trans2 = $ts->getTranscripts($genomic_pos2);
   
      $fusion = parse_transcript_data($fusion, 1, \@trans1);
      $fusion = parse_transcript_data($fusion, 2, \@trans2);
    }
    
		$breaklist{$fusion->{'breakpoint'}} = $fusion if(!exists $breaklist{$fusion->{'breakpoint'}});
  }
  close ($ifh2);
	
	my $final_annot_file1 = File::Spec->catfile($tmp, "$sample.1.ann_final");
	my $final_annot_file2 = File::Spec->catfile($tmp, "$sample.2.ann_final");
	my $bed_file1 = File::Spec->catfile($tmp, "$sample.1.bed");
	my $bed_file2 = File::Spec->catfile($tmp, "$sample.2.bed");
	my $final_annotation_file = File::Spec->catfile($tmp, "$sample.final");
	open(my $ofh1, '>', $final_annotation_file) or die "Could not open file '$final_annotation_file' $!";
	open(my $ofh2, '>', $final_annot_file1) or die "Could not open file '$final_annot_file1' $!";
	open(my $ofh3, '>', $final_annot_file2) or die "Could not open file '$final_annot_file2' $!";
	open(my $ofh4, '>', $bed_file1) or die "Could not open file '$bed_file1' $!";
	open(my $ofh5, '>', $bed_file2) or die "Could not open file '$bed_file2' $!";
	for my $brk (keys %breaklist){
	  if(defined $breaklist{$brk}->{'transcript1_id'} && defined $breaklist{$brk}->{'transcript2_id'}){
		  my $output_line = $breaklist{$brk}->format_annotation_line($vagrent_version);
		  print $ofh1 $output_line."\n";
		}
		else{
		  if(!defined $breaklist{$brk}->{'transcript1_id'}){
		    my $output_line = $breaklist{$brk}->format_bed_line(1);
		    print $ofh4 $output_line."\n";
		  }
		  else{
		    my $output_line = $breaklist{$brk}->format_break_line(1, $vagrent_version);
		    print $ofh2 $output_line."\n";
		  }
		  if(!defined $breaklist{$brk}->{'transcript2_id'}){
		    my $output_line = $breaklist{$brk}->format_bed_line(2);
		    print $ofh5 $output_line."\n";
		  }
		  else{
		    my $output_line = $breaklist{$brk}->format_break_line(2, $vagrent_version);
		    print $ofh3 $output_line."\n";
		  }
		}
  }
	close($ofh5);
	close($ofh4);
	close($ofh3);
	close($ofh2);
	close($ofh1);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
  return 1;
}

sub run_bed_pairtopair {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $sample = $options->{'sample'};
  
  my $prog = _which('bedtools');
  my @commands;
  
  # There will always be at least two input files so build the command for the first comparison
  my $overlap_file1 = File::Spec->catfile($tmp, "1_2.$sample.bedpe_overlap");
  push @commands, $prog . sprintf $BEDTOOLS_PAIRTOPAIR, 	File::Spec->catfile($tmp, "1.$sample.bedpe"),
  						 	File::Spec->catfile($tmp, "2.$sample.bedpe"),
  						 	$overlap_file1;
  						 	
  if($options->{'num'} == 3){
    
    my $overlap_file2 = File::Spec->catfile($tmp, "1_3.$sample.bedpe_overlap");
    push @commands, $prog . sprintf $BEDTOOLS_PAIRTOPAIR, 	File::Spec->catfile($tmp, "1.$sample.bedpe"),
  						 	File::Spec->catfile($tmp, "3.$sample.bedpe"),
  						 	$overlap_file2;
  	
  	my $overlap_file3 = File::Spec->catfile($tmp, "2_3.$sample.bedpe_overlap");				 	
  	push @commands, $prog . sprintf $BEDTOOLS_PAIRTOPAIR, 	File::Spec->catfile($tmp, "2.$sample.bedpe"),
  						 	File::Spec->catfile($tmp, "3.$sample.bedpe"),
  						 	$overlap_file3;
  						 	
  	 # Also find the intersection between all three algorithms using Unix sort and join commands.
  	 my $sorted_file1 = File::Spec->catfile($tmp, "1_2.$sample.bedpe_overlap.sorted");
  	 my $sorted_file2 = File::Spec->catfile($tmp, "1_3.$sample.bedpe_overlap.sorted");
  	 my $overlap_file123 = File::Spec->catfile($tmp, "1_2_3.$sample.bedpe_overlap");
  	 push @commands, sprintf $SORT, 8, 8, $overlap_file1, $sorted_file1;
  	 push @commands, sprintf $SORT, 8,8, $overlap_file2, $sorted_file2;
     push @commands, sprintf $JOIN, $sorted_file1,$sorted_file2, $overlap_file123;
  }
  
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

  return 1;
}

sub select_annotation {

  # All possible exon annotations have been retrieved for each breakpoint, we need to select annotation for the nearest.
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $sample = $options->{'sample'};
  
  my ($gtf, $path) = fileparse($options->{'gtf'});
	
  my $annot_file1 = File::Spec->catfile($tmp, "$sample.1.ann");
  my $annot_file2 = File::Spec->catfile($tmp, "$sample.2.ann");
  my $bed_file1 = File::Spec->catfile($tmp, "$sample.1.bed");
  my $bed_file2 = File::Spec->catfile($tmp, "$sample.2.bed");
  my $final_annot_file1 = $annot_file1."_final";
  my $final_annot_file2 = $annot_file2."_final";
  
  my $exon_annotation1;
  my $exon_annotation2;

	if(-s $annot_file1){
    $exon_annotation1 = process_annotation_file($annot_file1);
	}
	if(-s $annot_file2){
    $exon_annotation2 = process_annotation_file($annot_file2);
	}
	
	my $gene_gtf = File::Spec->catfile($tmp, "filtered_gene.gtf");
  my %gene_info;

  open (my $ifh1, $gene_gtf) or die "Could not open file '$gene_gtf' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    my $gene_annot = parse_gene_info($line);
    if(!exists $gene_info{$gene_annot->{'gene_name'}}){
      $gene_info{$gene_annot->{'gene_name'}}{'feature_start'} = $gene_annot->{'start'};
      $gene_info{$gene_annot->{'gene_name'}}{'feature_end'} = $gene_annot->{'end'};
    }
  }
  close ($ifh1);
	
	my $bed1;
  my $bed2;
  
  if(-s $bed_file1){
    open(my $ofh1, '>>', $final_annot_file1) or die "Could not open file '$final_annot_file1' $!";
    open (my $ifh2, $bed_file1) or die "Could not open file '$bed_file1' $!";
    while (<$ifh2>) {
      chomp;
      my $line = $_;
      my $fusion = parse_bed_file($line, 1);
      if(exists $exon_annotation1->{$fusion->breakpoint}){
        $fusion = parse_exon_data($fusion, 1, $exon_annotation1->{$fusion->breakpoint});
      }
      if(!exists $fusion->{'transcript1_id'}){
        $fusion->{'exon1_num'} = 'NA';
        $fusion->{'feature1_start'} = 'NA';
        $fusion->{'feature1_end'} = 'NA';
        # Check whether the breakpoint falls within the footprint of the gene.
        my $gene_start = $gene_info{$fusion->{'gene1'}}{'feature_start'};
			  my $gene_end = $gene_info{$fusion->{'gene1'}}{'feature_end'};
			  my $break_pos = $fusion->{'pos1_end'};
			  if(defined $gene_start && defined $gene_end){
			    if($break_pos >= $gene_start && $break_pos <= $gene_end){
			      $fusion->{'transcript1_id'} = 'Intronic';
			    }
			    else{
			      my $location = 'upstream';
			      $location = 'downstream' if($break_pos > $gene_end);
			  		$fusion->{'transcript1_id'} = $location;
			    }
			  }
			  else{
			    $fusion->{'transcript1_id'} = 'Unannotated';
			  }
      }
      my $output_line = $fusion->format_break_line(1, $gtf);
		  print $ofh1 $output_line."\n";
    }
    close($ifh2);
    close($ofh1);
	}
	
	if(-s $bed_file2){
	  open(my $ofh2, '>>', $final_annot_file2) or die "Could not open file '$final_annot_file2' $!";
    open (my $ifh3, $bed_file2) or die "Could not open file '$bed_file2' $!";
    while (<$ifh3>) {
      chomp;
      my $line = $_;
      my $fusion = parse_bed_file($line, 2);
      if(exists $exon_annotation1->{$fusion->breakpoint}){
        $fusion = parse_exon_data($fusion, 2, $exon_annotation1->{$fusion->breakpoint});
      }
      if(!exists $fusion->{'transcript2_id'}){
        $fusion->{'exon2_num'} = 'NA';
        $fusion->{'feature2_start'} = 'NA';
        $fusion->{'feature2_end'} = 'NA';
        # Check whether the breakpoint falls within the footprint of the gene.
        my $gene_start = $gene_info{$fusion->{'gene2'}}{'feature_start'};
			  my $gene_end = $gene_info{$fusion->{'gene2'}}{'feature_end'};
			  my $break_pos = $fusion->{'pos2_end'};
			  if(defined $gene_start && defined $gene_end){
			    if($break_pos >= $gene_start && $break_pos <= $gene_end){
			      $fusion->{'transcript2_id'} = 'Intronic';
			    }
			    else{
			      my $location = 'upstream';
			      $location = 'downstream' if($break_pos > $gene_end);
			  		$fusion->{'transcript2_id'} = $location;
			    }
			  }
			  else{
			    $fusion->{'transcript2_id'} = 'Unannotated';
			  }
      }
      my $output_line = $fusion->format_break_line(2, $gtf);
		  print $ofh2 $output_line."\n";
    }
    close($ifh3);
    close($ofh2);
	}

  #PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
  return 1;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in path or local bin folder ($l_bin)\n\tPATH: $ENV{PATH}\n" unless(defined $path && -e $path);
  return $path;
}

1;
__END__

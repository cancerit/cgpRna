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

const my $OUTPUT_HEADER => "sample\tfusion_name\tstar_junction\tdefuse_junction\tdefuse_splitr_count\tdefuse_span_count\tstar_JunctionReads\tstar_SpanningFrags\tLeftGene\tLeftGeneId\tLeftChr\tLeftPos\tLeftStrand\tLeftDistFromRefExonSplice\tRightGene\tRightGeneId\tRightChr\tRightPos\tRightStrand\tRightDistFromRefExonSplice\ttranscript1_id\ttranscript1_src\texon1_num\texon1_start\texon1_end\ttranscript2_id\ttranscript2_src\texon2_num\texon2_start\texon2_end\tdefuse_cluster_id\tdefuse_splitr_sequence\n";

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

# Position of the columns in the tophat-fusion filtered file used to format the bed file.
const my $TOPHAT_SPLIT_CHAR => '\t';
const my $TOPHAT_CHR1 => 4;
const my $TOPHAT_POS1 => 5;
const my $TOPHAT_STRAND1 => 13;
const my $TOPHAT_CHR2 => 7;
const my $TOPHAT_POS2 => 8;
const my $TOPHAT_STRAND2 => 14;
const my $TOPHAT_BREAKREF => 1;
const my $TOPHAT_HEADER_PATTERN => 'num_spanning_reads';

# Position of the columns in the deFuse output file used to format fusion breakpoint references.
const my $DEFUSE_SPLIT_CHAR => '\t';
const my $DEFUSE_CHR1 => 26;
const my $DEFUSE_POS1 => 39;
const my $DEFUSE_STRAND1 => 36;
const my $DEFUSE_CHR2 => 27;
const my $DEFUSE_POS2 => 40;
const my $DEFUSE_STRAND2 => 37;
const my $DEFUSE_BREAKREF => 1;
const my $DEFUSE_CLUSTER_ID => 2;
const my $DEFUSE_HEADER_PATTERN => 'cluster_id';

# Position of the columns in the star-fusion output file used to format fusion breakpoint references.
const my $STAR_SPLIT_CHAR => '\t';
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
const my $STAR_HEADER_PATTERN => 'fusion_name';

sub annotate_bed {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $sample = $options->{'sample'};
  my $exon_gtf = filter_gtf($options, 'exon');
  my $gene_gtf = filter_gtf($options, 'gene');

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
    my $breakpoint_id = $fields[0];
    open (my $ifh2, $annot_file2) or die "Could not open file '$annot_file2' $!";
    while (<$ifh2>) {
      chomp;
      my $line2 = $_;
      if($line2 =~ m/^$breakpoint_id/){
        my @fields2 = split "\t", $line2;
        shift(@fields2);
				shift(@fields2);
				push(@fields, @fields2);
				last;
      }
    }
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
    print $ofh $line."\n" if(exists $ALLOWED_BIOTYPES{$annotation{'gene_biotype'}});
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
  
  my $star_file;
  my $defuse_file;
  
  if($options->{'fusion_files'}->{'1'}->{'format'} eq 'star'){
    $star_file = $options->{'fusion_files'}->{'1'}->{'name'};
    $defuse_file = $options->{'fusion_files'}->{'2'}->{'name'};
  }
  else{
    $star_file = $options->{'fusion_files'}->{'2'}->{'name'};
    $defuse_file = $options->{'fusion_files'}->{'1'}->{'name'};
  }

  my %star_data;
  open (my $ifh1, $star_file) or die "Could not open file '$star_file' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    next if($line =~ m/^breakpoint/);
    $line =~ m/^(.*:[0-9]+-.*:[0-9]+)\t([A-Za-z0-9-_:\.]+--[A-Za-z0-9-_:\.]+)\t/;
    my $break_ref = $1;
		my $fusion_name = $2;
		$line =~ s/^.*:[0-9]+-.*:[0-9]+\t[A-Za-z0-9-_:\.]+--[A-Za-z0-9-_:\.]+\t//;
		$star_data{$break_ref}{'fusion_name'} = $fusion_name;
    $star_data{$break_ref}{'data'} = $line;
  }
  close ($ifh1);
  
  my %defuse_data;
  open (my $ifh2, $defuse_file) or die "Could not open file '$defuse_file' $!";
  while (<$ifh2>) {
    chomp;
    my $line = $_;
    next if($line =~ m/^breakpoint/);
    my @fields = split "\t", $line;
    $line =~ m/^(.*:[0-9]+-.*:[0-9]+)\t([0-9]+)\t([ACGT|]+)\t/;
    my $break_ref = $fields[0];
    my $cluster_id = $fields[1];
    my $sequence = $fields[2];
    my $split_reads = $fields[3];
    my $span_reads = $fields[61];
    $defuse_data{$break_ref."_".$cluster_id}{'breakpoint'} = $break_ref;
    $defuse_data{$break_ref."_".$cluster_id}{'cluster_id'} = $cluster_id;
    $defuse_data{$break_ref."_".$cluster_id}{'sequence'} = $sequence;
    $defuse_data{$break_ref."_".$cluster_id}{'split_reads'} = $split_reads;
    $defuse_data{$break_ref."_".$cluster_id}{'span_reads'} = $span_reads;
  }
  close ($ifh2);

  my $annot_file = File::Spec->catfile($tmp, "$sample.final");
  my $output_file = File::Spec->catfile($tmp, "$sample.star-defuse.overlapping.fusions.txt");
  open(my $ofh1, '>', $output_file) or die "Could not open file $output_file $!";
  
  if(-s $annot_file){
  open(my $ifh3, $annot_file) or die "Could not open file $annot_file $!";
  print $ofh1 $OUTPUT_HEADER;
  while(<$ifh3>){
    chomp;
    my $line = $_;
    my @fields = split "\t", $line;
    my $breakpoint = $fields[0];
    my $alt_breakpoint = $fields[1];
    my $defuse_breakpoint = $defuse_data{$alt_breakpoint}{'breakpoint'};
  	my $length = scalar @fields;
  	if($length > 10){
      my $star_fusion_name = $star_data{$breakpoint}{'fusion_name'};
      my $star_data = $star_data{$breakpoint}{'data'};
      my $exon1_number = $fields[6];
      my $exon2_number = $fields[13];
      my $exon1_start = $fields[7];
      my $exon2_start = $fields[14];
      my $exon1_end = $fields[8];
      my $exon2_end = $fields[15];
      my $transcript1_id = $fields[4];
      my $transcript2_id = $fields[11];
      my $transcript1_src = $fields[5];
      my $transcript2_src = $fields[12];
      my $defuse_cluster_id = $defuse_data{$alt_breakpoint}{'cluster_id'};
      my $defuse_split_reads = $defuse_data{$alt_breakpoint}{'split_reads'};
      my $defuse_span_reads = $defuse_data{$alt_breakpoint}{'span_reads'};
      my $defuse_sequence = $defuse_data{$alt_breakpoint}{'sequence'};
      print $ofh1 "$sample\t$star_fusion_name\t$breakpoint\t$defuse_breakpoint\t$defuse_split_reads\t$defuse_span_reads\t$star_data\t$transcript1_id\t$transcript1_src\t$exon1_number\t$exon1_start\t$exon1_end\t$transcript2_id\t$transcript2_src\t$exon2_number\t$exon2_start\t$exon2_end\t$defuse_cluster_id\t$defuse_sequence\n";
    }
    else{
      print $ofh1 "$sample\t-\t$breakpoint\t$defuse_breakpoint\tFUSION COULD NOT BE ANNOTATED\n";
    }
  }
  close($ifh3);
  }
  else{
    print $ofh1 "$sample\tNO OVERLAPPING FUSIONS FOUND\n";
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
  $annotation{'chr'} = $fields[0];
  $annotation{'pos_start'} = $fields[1];
  $annotation{'pos_end'} = $fields[2];
  $annotation{'strand'} = $fields[5];
  $annotation{'feature'} = $fields[10];
  $annotation{'feature_start'} = $fields[11];
  $annotation{'feature_end'} = $fields[12];
  $annotation{'star_genename'} = $fields[6];
  $annotation{'star_geneid'} = $fields[7];

	my $annot_column = scalar @fields;
  my @annot_fields = split /; /, $fields[$annot_column-1];
  foreach my $item(@annot_fields) {
    my ($type,$value)= split / /, $item;
    $annotation{$type} = $value;
  }

  return \%annotation;
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

sub parse_overlap {
  my ($line, $cols) = @_;
  
  my @fields = split "\t", $line;
  
  my $row_length = scalar @fields;
  my $start = 0;
  $start = $row_length / 2  if($cols == 2);
  
  my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
    -breakpoint	=> $fields[$start + 7],
    -chr1	=> $fields[$start],
    -pos1_start	=> $fields[$start + 1],
    -pos1_end	=> $fields[$start + 2],
    -strand1	=> $fields[$start + 8],
    -chr2	=> $fields[$start + 3],
    -pos2_start	=> $fields[$start + 4],
    -pos2_end	=> $fields[$start + 5],
    -strand2	=> $fields[$start + 9]);
   
  if($cols == 1){
    $fusion->alt_breakpoint($fields[17]);
  }
  else{
    $fusion->alt_breakpoint($fields[7]);
  }
  
  return $fusion;
}

sub parse_transcript_data {
  my ($fusion, $breaknum, $transcripts) = @_;
  
  my @filteredTrans;
  foreach my $t(@{$transcripts}){
		push(@filteredTrans, $t) if ($fusion->{'pos'.$breaknum.'_end'} >= $t->getGenomicMinPos && $fusion->{'pos'.$breaknum.'_end'} <= $t->getGenomicMaxPos);
	}
	my @sortedTrans = sort{&annotation_sort} @filteredTrans;
	
	if(defined $sortedTrans[0]){

	  my $exon_number;
		my $exon_start;
		my $exon_end;
		my $transcript_id;
		my $gene_biotype;
		my $num_transcripts = scalar @sortedTrans;
		for (my $x=0;$x<$num_transcripts; $x++){
		  my @exons = $sortedTrans[$x]->getExons;
		  my $num_exons = scalar @exons;
		  for (my $y=0;$y<$num_exons; $y++){
		    my $e = $exons[$y];
		    if($fusion->{'pos'.$breaknum.'_end'} == $e->getMinPos || $fusion->{'pos'.$breaknum.'_end'} == $e->getMaxPos){
		      $transcript_id = $sortedTrans[$x]->getAccession;
		      $gene_biotype = $sortedTrans[$x]->{'_genetype'};
		      $exon_start = $e->getMinPos;
		      $exon_end = $e->getMaxPos;
		      $exon_number = $y+1;
		      last;
		    }
		  }
		  last if(defined $transcript_id);
		}
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

sub process_annotation_file {
  my ($options, $input, $output, $gene_info) = @_;

  my $curr_distance = 10000000;
  my $curr_break = "";
  my $curr_annotation;
  my $curr_pos;
  my $curr_exon_start;
	my $curr_exon_end;
  
  open(my $ofh1, '>>', $output) or die "Could not open file $output $!";
  open (my $ifh1, $input) or die "Could not open file '$input' $!";
  while (<$ifh1>){
    chomp;
    my $line = $_;
    my $annotation = parse_annotation($line);
		next if($annotation->{'gene_name'} ne $annotation->{'star_genename'});
		my $break = $annotation->{'breakpoint'};
		if($break ne $curr_break){
		  unless($curr_break eq ""){
		    $curr_pos = $curr_annotation->{'pos_end'};
		    $curr_exon_start = $curr_annotation->{'feature_start'};
    		$curr_exon_end = $curr_annotation->{'feature_end'};
    		
		    if($curr_distance <= 10){
          print $ofh1 $curr_annotation->{'breakpoint'}."\t".$curr_annotation->{'alt_breakpoint'}."\t".$curr_annotation->{'star_genename'}."\t".$curr_annotation->{'star_geneid'}."\t".$curr_annotation->{'transcript_id'}."\t".$options->{'gtf'}."\t".$curr_annotation->{'exon_number'}."\t".$curr_annotation->{'feature_start'}."\t".$curr_annotation->{'feature_end'}."\n";
			  }
			  elsif($curr_pos > $curr_exon_start && $curr_pos < $curr_exon_end){
          print $ofh1 $curr_annotation->{'breakpoint'}."\t".$curr_annotation->{'alt_breakpoint'}."\t".$curr_annotation->{'star_genename'}."\t".$curr_annotation->{'star_geneid'}."\t".$curr_annotation->{'transcript_id'}."\t".$options->{'gtf'}."\t".$curr_annotation->{'exon_number'}."\t".$curr_annotation->{'feature_start'}."\t".$curr_annotation->{'feature_end'}."\n";			    
			  }
			  else{
			    # The breakpoint doesn't fall within 10bp of an exon boundary. We need to check it falls within the footprint of the star gene and, for now, print it as intronic
			    my $gene_start = $gene_info->{$curr_annotation->{'star_genename'}}{'feature_start'};
			    my $gene_end = $gene_info->{$curr_annotation->{'star_genename'}}{'feature_end'};
			    my $break_pos = $curr_annotation->{'pos_end'};
			    if($break_pos >= $gene_start && $break_pos <= $gene_end){
			      print $ofh1 $curr_annotation->{'breakpoint'}."\t".$curr_annotation->{'alt_breakpoint'}."\t".$curr_annotation->{'star_genename'}."\t".$curr_annotation->{'star_geneid'}."\tIntronic\t".$options->{'gtf'}."\t-\t-\t-\n";
			    }
			  }
			  $curr_distance = 10000000;
			  $curr_break = $break;
		  }
		}
		my $pos = $annotation->{'pos_end'};
    my $exon_start = $annotation->{'feature_start'};
    my $exon_end = $annotation->{'feature_end'};
    my $distance = find_closest_boundary($pos, $exon_start, $exon_end);

    if($distance < $curr_distance){
      $curr_distance = $distance;
      $curr_annotation = $annotation;
      $curr_break = $break;
    }
  }
  $curr_pos = $curr_annotation->{'pos_end'};
	$curr_exon_start = $curr_annotation->{'feature_start'};
  $curr_exon_end = $curr_annotation->{'feature_end'};
  
	if($curr_distance <= 10){
   print $ofh1 $curr_annotation->{'breakpoint'}."\t".$curr_annotation->{'alt_breakpoint'}."\t".$curr_annotation->{'star_genename'}."\t".$curr_annotation->{'star_geneid'}."\t".$curr_annotation->{'transcript_id'}."\t".$options->{'gtf'}."\t".$curr_annotation->{'exon_number'}."\t".$curr_annotation->{'feature_start'}."\t".$curr_annotation->{'feature_end'}."\n";
	}
  elsif($curr_pos > $curr_exon_start && $curr_pos < $curr_exon_end){
    print $ofh1 $curr_annotation->{'breakpoint'}."\t".$curr_annotation->{'alt_breakpoint'}."\t".$curr_annotation->{'star_genename'}."\t".$curr_annotation->{'star_geneid'}."\t".$curr_annotation->{'transcript_id'}."\t".$options->{'gtf'}."\t".$curr_annotation->{'exon_number'}."\t".$curr_annotation->{'feature_start'}."\t".$curr_annotation->{'feature_end'}."\n";			    
  }
	else{
	  # The breakpoint doesn't fall within 10bp of an exon boundary. We need to check it falls within the footprint of the star gene and, for now, print it as intronic
	  my $gene_start = $gene_info->{$curr_annotation->{'star_genename'}}{'feature_start'};
	  my $gene_end = $gene_info->{$curr_annotation->{'star_genename'}}{'feature_end'};
	  my $break_pos = $curr_annotation->{'pos_end'};
	  if($break_pos >= $gene_start && $break_pos <= $gene_end){
	    print $ofh1 $curr_annotation->{'breakpoint'}."\t".$curr_annotation->{'alt_breakpoint'}."\t".$curr_annotation->{'star_genename'}."\t".$curr_annotation->{'star_geneid'}."\tIntronic\t".$options->{'gtf'}."\t-\t-\t-\n";
		}
  }
  close ($ifh1);
  close ($ofh1);

  return 1;
}

sub query_vagrent {
  my $options = shift;

  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
  my $sample = $options->{'sample'};
  
  # There will always be a 1_2 comparison file so deal with that first and build the fusions object.
  
  # Establish the source of 1 and 2 respectively
  my $source1 = $options->{'fusion_files'}->{'1'}->{'format'};
  my $source2 = $options->{'fusion_files'}->{'2'}->{'format'};
  
  my $col_set = 1;
  my $star_file = $options->{'fusion_files'}->{'1'}->{'name'};
  if($source2 eq 'star'){
    $col_set = 2;
    $star_file = $options->{'fusion_files'}->{'2'}->{'name'};
  }

  my $overlap_file1_2 = File::Spec->catfile($tmp, "1_2.$sample.bedpe_overlap");
  
  my %star_gene_list;
  open (my $ifh1, $star_file) or die "Could not open file '$star_file' $!";
  while (<$ifh1>) {
    chomp;
    my $line = $_;
    my @fields = split "\t", $line;
    my $breakpoint = $fields[0];
    $star_gene_list{$breakpoint}{'gene1_name'} = $fields[4];
    $star_gene_list{$breakpoint}{'gene1_id'} = $fields[5];
    $star_gene_list{$breakpoint}{'gene2_name'} = $fields[10];
    $star_gene_list{$breakpoint}{'gene2_id'} = $fields[11];
  }
  close ($ifh1);
  
	my $vagrent_version = "VAGrENT_".Sanger::CGP::Vagrent->VERSION;
	
	my $ts = Sanger::CGP::Vagrent::TranscriptSource::FileBasedTranscriptSource->new('cache' => $options->{'cache'});	
	my %breaklist;
	
  open (my $ifh3, $overlap_file1_2) or die "Could not open file '$overlap_file1_2' $!";
  while (<$ifh3>) {
    chomp;
    my $line = $_;
    my $fusion = parse_overlap($line, $col_set);
    $fusion->gene1($star_gene_list{$fusion->{'breakpoint'}}{'gene1_name'});
    $fusion->gene2($star_gene_list{$fusion->{'breakpoint'}}{'gene2_name'});
    $fusion->gene1_id($star_gene_list{$fusion->{'breakpoint'}}{'gene1_id'});
    $fusion->gene2_id($star_gene_list{$fusion->{'breakpoint'}}{'gene2_id'});
    
    my $genomic_pos1 = Sanger::CGP::Vagrent::Data::GenomicRegion->new('species' => 'human', 'genomeVersion' => 'GRCh38', 'chr' => $fusion->{'chr1'}, 'minpos' => $fusion->{'pos1_start'}, 'maxpos' => $fusion->{'pos1_end'}, 'id' => $fusion->{'breakpoint'});
    my $genomic_pos2 = Sanger::CGP::Vagrent::Data::GenomicRegion->new('species' => 'human', 'genomeVersion' => 'GRCh38', 'chr' => $fusion->{'chr2'}, 'minpos' => $fusion->{'pos2_start'}, 'maxpos' => $fusion->{'pos2_end'}, 'id' => $fusion->{'breakpoint'});
		
    my @trans1 = $ts->getTranscripts($genomic_pos1);
    my @trans2 = $ts->getTranscripts($genomic_pos2);
    
    parse_transcript_data($fusion, 1, \@trans1);
    parse_transcript_data($fusion, 2, \@trans2);

		$breaklist{$fusion->{'breakpoint'}} = $fusion if(!exists $breaklist{$fusion->{'breakpoint'}});
  }
  close ($ifh3);
	
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
  
  # There will always be at least two input files so build the command for the first comparison
  my $command = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, 	File::Spec->catfile($tmp, "1.$sample.bedpe"),
  						 	File::Spec->catfile($tmp, "2.$sample.bedpe"),
  						 	File::Spec->catfile($tmp, "1_2.$sample.bedpe_overlap");
  
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

  return 1;
}

sub select_annotation {

  # All possible exon annotations have been retrieved for each breakpoint, we need to select annotation for the nearest.
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $sample = $options->{'sample'};
  my $annot_file1;
  my $annot_file2;
  my $final_annot_file1;
  my $final_annot_file2;
  
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
	
  opendir(my $dh, $tmp);
  while(my $file = readdir $dh){
    $annot_file1 = File::Spec->catfile($tmp, $file) if($file eq "$sample.1.ann");
    $annot_file2 = File::Spec->catfile($tmp, $file) if($file eq "$sample.2.ann");
  }
  closedir($dh);
	
  $final_annot_file1 = $annot_file1."_final";
  $final_annot_file2 = $annot_file2."_final";

	if(-s $annot_file1){
    process_annotation_file($options, $annot_file1, $final_annot_file1, \%gene_info);
	}
	if(-s $annot_file2){
    process_annotation_file($options, $annot_file2, $final_annot_file2, \%gene_info);
	}	
    
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
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

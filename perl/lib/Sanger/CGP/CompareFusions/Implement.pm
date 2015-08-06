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

use Data::Dumper;

const my $BEDTOOLS_CLOSEST => q{ closest -d -a %s -b %s -t first > %s};
const my $BEDTOOLS_PAIRTOPAIR => q{ pairtopair -a %s -b %s -f 1.0 > %s};

const my $OUTPUT_GENE_HEADER => "genes\tbreakpoints\tcalled_by\tchr1\tstrand1\tgene1\tgene1_start\tgene1_end\tchr2\tstrand2\tgene2\tgene2_start\tgene2_end\n";
const my $OUTPUT_EXON_HEADER => "exons\tbreakpoints\tcalled_by\tchr1\tstrand1\tgene1\ttranscript1_id\texon1_num\texon1_start\texon1_end\tchr2\tstrand2\tgene2\ttranscript2_id\texon2_num\texon2_start\texon2_end\n";

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
const my $STAR_HEADER_PATTERN => 'fusion_name';

# Position of the columns in the SOAPfuse output file used to format fusion breakpoint references.
const my $SOAP_SPLIT_CHAR => '\t';
const my $SOAP_CHR1 => 3;
const my $SOAP_POS1 => 5;
const my $SOAP_STRAND1 => 4;
const my $SOAP_CHR2 => 8;
const my $SOAP_POS2 => 10;
const my $SOAP_STRAND2 => 9;
const my $SOAP_BREAKREF => 1;
const my $SOAP_HEADER_PATTERN => 'up_chr';

sub annotate_bed {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});
  
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $sample = $options->{'sample'};
	my $exon_gtf = $options->{'exon_gtf'};
	my $gene_gtf = $options->{'gene_gtf'};

	my $break1_file;
	my $break2_file;
	
	opendir(my $dh, $tmp);
	while(my $file = readdir $dh) {
	  $break1_file = File::Spec->catfile($tmp, $file) if($file =~ m/^$index.*1.bed/);
	  $break2_file = File::Spec->catfile($tmp, $file) if($file =~ m/^$index.*2.bed/);
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

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, $index);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);

  return 1;
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
	# SOAPfuse check
	elsif($firstLine =~ m/$SOAP_HEADER_PATTERN/){
		$source = "soap";
	}
	else{
		die "Unrecognised file type or the file is missing the header record\n";
	}
	
	return $source;
}

sub compare_overlaps {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $sample = $options->{'sample'};
  
  # There will always be a 1_2 comparison file so deal with that first and build the fusions object.
  
  # Establish the source of 1 and 2 respectively
  my $source1 = $options->{'fusion_files'}->{'1'}->{'format'};
  my $source2 = $options->{'fusion_files'}->{'2'}->{'format'};
  
  my $gene_overlap_file1_2;
  my $exon_overlap_file1_2;
  
  opendir(my $dh, $tmp);
	while(my $file = readdir $dh) {
	  $gene_overlap_file1_2 = File::Spec->catfile($tmp, $file) if($file =~ m/^1_2.$sample.gene.bedpe_overlap/);
	  $exon_overlap_file1_2 = File::Spec->catfile($tmp, $file) if($file =~ m/^1_2.$sample.exon.bedpe_overlap/);
	}
	closedir($dh);
	
	my %gene_list;
	my %feature_source;
	my %exon_list;
	
  process_gene_overlaps($gene_overlap_file1_2, \%gene_list, \%feature_source, $source1, $source2);
	process_exon_overlaps($exon_overlap_file1_2, \%exon_list, \%feature_source, $source1, $source2);
		
  if($options->{'num'} == 3){
    my $gene_overlap_file1_3;
    my $exon_overlap_file1_3;
    my $gene_overlap_file2_3;
    my $exon_overlap_file2_3;
    
    my $source3 = $options->{'fusion_files'}->{'3'}->{'format'};
    
    opendir(my $dh2, $tmp);
	  while(my $file = readdir $dh2) {
	    $gene_overlap_file1_3 = File::Spec->catfile($tmp, $file) if($file =~ m/^1_3.$sample.gene.bedpe_overlap/);
	    $exon_overlap_file1_3 = File::Spec->catfile($tmp, $file) if($file =~ m/^1_3.$sample.exon.bedpe_overlap/);
	    $gene_overlap_file2_3 = File::Spec->catfile($tmp, $file) if($file =~ m/^2_3.$sample.gene.bedpe_overlap/);
	    $exon_overlap_file2_3 = File::Spec->catfile($tmp, $file) if($file =~ m/^2_3.$sample.exon.bedpe_overlap/);
	  }
	  closedir($dh2);
	  
	  process_gene_overlaps($gene_overlap_file1_3, \%gene_list, \%feature_source, $source1, $source3);
	  process_exon_overlaps($exon_overlap_file1_3, \%exon_list, \%feature_source, $source1, $source3);
	  process_gene_overlaps($gene_overlap_file2_3, \%gene_list, \%feature_source, $source2, $source3);
	  process_exon_overlaps($exon_overlap_file2_3, \%exon_list, \%feature_source, $source2, $source3);
  
  }

  my $gene_output_file = File::Spec->catfile($tmp, "$sample.gene-fusions.txt");
  open(my $ofh1, '>', $gene_output_file) or die "Could not open file $gene_output_file $!";
  print $ofh1 $OUTPUT_GENE_HEADER;
  for my $gene (keys %gene_list){
    my @brk_list;
    for my $brk (keys %{ $gene_list{$gene}}){
      push @brk_list, $brk;
  	}
  	
  	my @brks;
  	for my $brkl (keys %{ $feature_source{$gene}{'breakpoints'}}){
      push @brks, $brkl;
  	}

  	my $tophat = "";
  	my $defuse = "";
  	my $star = "";
  	my $soap = "";
  	
  	$tophat = "T" if(exists $feature_source{$gene}->{'tophat'});
  	$defuse = "D" if(exists $feature_source{$gene}->{'defuse'});
  	$star = "S" if(exists $feature_source{$gene}->{'star'});
    $soap = "O" if(exists $feature_source{$gene}->{'soap'});
  	my $gene_source_string = $tophat.$defuse.$star.$soap;  	
  	
  	my $chr1 = $gene_list{$gene}{$brk_list[0]}->{'chr1'};
  	my $strand1 = $gene_list{$gene}{$brk_list[0]}->{'strand1'};
  	my $gene1 = $gene_list{$gene}{$brk_list[0]}->{'gene1'};
  	my $gene1_start = $gene_list{$gene}{$brk_list[0]}->{'gene1_start'};
  	my $gene1_end = $gene_list{$gene}{$brk_list[0]}->{'gene1_end'};
  	my $chr2 = $gene_list{$gene}{$brk_list[0]}->{'chr2'};
  	my $strand2 = $gene_list{$gene}{$brk_list[0]}->{'strand2'};
  	my $gene2 = $gene_list{$gene}{$brk_list[0]}->{'gene2'};
  	my $gene2_start = $gene_list{$gene}{$brk_list[0]}->{'gene2_start'};
  	my $gene2_end = $gene_list{$gene}{$brk_list[0]}->{'gene2_end'};
  	  
    my $breaks = join(",",@brks);
    print $ofh1 "$gene\t$breaks\t$gene_source_string\t$chr1\t$strand1\t$gene1\t$gene1_start\t$gene1_end\t$chr2\t$strand2\t$gene2\t$gene2_start\t$gene2_end\n";
  }
  close($ofh1);
    
  my $exon_output_file = File::Spec->catfile($tmp, "$sample.exon-fusions.txt");
  open(my $ofh2, '>', $exon_output_file) or die "Could not open file $exon_output_file $!";
  print $ofh2 $OUTPUT_EXON_HEADER;
    
  for my $exon (keys %exon_list){
    my @brk_list;
    for my $brk (keys %{ $exon_list{$exon}}){
      push @brk_list, $brk;
    }
    
    my @brks;
  	for my $brkl (keys %{ $feature_source{$exon}{'breakpoints'}}){
      push @brks, $brkl;
  	}

  	my $tophat = "";
  	my $defuse = "";
  	my $star = "";
  	my $soap = "";
    $tophat = "T" if(exists $feature_source{$exon}->{'tophat'});
    $defuse = "D" if(exists $feature_source{$exon}->{'defuse'});
    $star = "S" if(exists $feature_source{$exon}->{'star'});
    $soap = "O" if(exists $feature_source{$exon}->{'soap'});
    my $exon_source_string = $tophat.$defuse.$star.$soap;   
    
    my $chr1 = $exon_list{$exon}{$brk_list[0]}->{'chr1'};
    my $strand1 = $exon_list{$exon}{$brk_list[0]}->{'strand1'};
    my $gene1 = $exon_list{$exon}{$brk_list[0]}->{'gene1'};
    my $transcript1_id = $exon_list{$exon}{$brk_list[0]}->{'transcript1_id'};
    my $exon1_num = $exon_list{$exon}{$brk_list[0]}->{'exon1_num'};
    my $exon1_start = $exon_list{$exon}{$brk_list[0]}->{'feature1_start'};
    my $exon1_end = $exon_list{$exon}{$brk_list[0]}->{'feature1_end'};
    my $chr2 = $exon_list{$exon}{$brk_list[0]}->{'chr2'};
    my $strand2 = $exon_list{$exon}{$brk_list[0]}->{'strand2'};
    my $gene2 = $exon_list{$exon}{$brk_list[0]}->{'gene2'};
    my $transcript2_id = $exon_list{$exon}{$brk_list[0]}->{'transcript2_id'};
    my $exon2_num = $exon_list{$exon}{$brk_list[0]}->{'exon2_num'};
    my $exon2_start = $exon_list{$exon}{$brk_list[0]}->{'feature2_start'};
    my $exon2_end = $exon_list{$exon}{$brk_list[0]}->{'feature2_end'};
  	  
    my $breaks = join(",",@brks);
    print $ofh2 "$exon\t$breaks\t$exon_source_string\t$chr1\t$strand1\t$gene1\t$transcript1_id\t$exon1_num\t$exon1_start\t$exon1_end\t$chr2\t$strand2\t$gene2\t$transcript2_id\t$exon2_num\t$exon2_start\t$exon2_end\n";
  }
  close($ofh2);

  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

  return 1;
}

sub create_bed {
	my ($index, $options) = @_;

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $sample = $options->{'sample'};
	
	my $file = $options->{'fusion_files'}->{$index}->{'name'};
	my $filetype = $options->{'fusion_files'}->{$index}->{'format'};
  my $output1 = File::Spec->catfile($tmp, "$index.$sample.$filetype.1.bed");
	my $output2 = File::Spec->catfile($tmp, "$index.$sample.$filetype.2.bed");
	
	open (my $ifh, $file) or die "Could not open file '$file' $!";
	open(my $ofh1, '>', $output1) or die "Could not open file '$output1' $!";
	open(my $ofh2, '>', $output2) or die "Could not open file '$output2' $!";
		
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
		  
		if($filetype eq 'tophat'){
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
		    		    
		  print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$name."\t".$strand1."\n";
		  print $ofh2 $chr2."\t".$pos2_start."\t".$pos2_end."\t".$name."\t".$strand2."\n";
		}
		elsif($filetype eq 'star'){
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
		    		    
		  print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$name."\t".$strand1."\n";
		  print $ofh2 $chr2."\t".$pos2_start."\t".$pos2_end."\t".$name."\t".$strand2."\n";
		  	
		}
		elsif($filetype eq 'soap'){
			next if($line =~ m/$SOAP_HEADER_PATTERN/);
		  	
		 	@fields = split $SOAP_SPLIT_CHAR, $line;
		  $name = $fields[$SOAP_BREAKREF - 1];
		  $chr1 = $fields[$SOAP_CHR1 - 1];
		  $pos1_start = $fields[$SOAP_POS1 - 1]-1;
		  $pos1_end = $fields[$SOAP_POS1 - 1];
		  $strand1 = $fields[$SOAP_STRAND1 - 1];
		  $chr2 = $fields[$SOAP_CHR2 - 1];
		  $pos2_start = $fields[$SOAP_POS2 - 1]-1;
		  $pos2_end = $fields[$SOAP_POS2- 1];
		  $strand2 = $fields[$SOAP_STRAND2 - 1];
		    		    
		  print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$name."\t".$strand1."\n";
		  print $ofh2 $chr2."\t".$pos2_start."\t".$pos2_end."\t".$name."\t".$strand2."\n";
		  	
		}
		# It must be defuse format
		else{
			next if($line =~ m/$DEFUSE_HEADER_PATTERN/);
		  	
			@fields = split $DEFUSE_SPLIT_CHAR, $line;
		  $name = $fields[$DEFUSE_BREAKREF - 1];
		  $chr1 = $fields[$DEFUSE_CHR1 - 1];
		  $pos1_start = $fields[$DEFUSE_POS1 - 1]-1;
		  $pos1_end = $fields[$DEFUSE_POS1 - 1];
		  $strand1 = $fields[$DEFUSE_STRAND1 - 1];
		  $chr2 = $fields[$DEFUSE_CHR2 - 1];
		  $pos2_start = $fields[$DEFUSE_POS2 - 1]-1;
		  $pos2_end = $fields[$DEFUSE_POS2- 1];
		  $strand2 = $fields[$DEFUSE_STRAND2 - 1];
		    		    
		  print $ofh1 $chr1."\t".$pos1_start."\t".$pos1_end."\t".$name."\t".$strand1."\n";
		  print $ofh2 $chr2."\t".$pos2_start."\t".$pos2_end."\t".$name."\t".$strand2."\n";
		}
		  
  }
	close ($ifh);
	close ($ofh1);
	close ($ofh2);
	
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
	
	return 1;
}

sub create_bedpe {
  my ($index, $options) = @_;
  return 1 if(exists $options->{'index'} && $index != $options->{'index'});

  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
  
	my $sample = $options->{'sample'};
	my $gtffile = $options->{'gtf'};
	
	my $annot_file1;
	my $annot_file2;
	my $gene_bedpe_file = File::Spec->catfile($tmp, "$index.$sample.gene.bedpe");
	my $exon_bedpe_file = File::Spec->catfile($tmp, "$index.$sample.exon.bedpe");
	
	opendir(my $dh, $tmp);
	while(my $file = readdir $dh) {
	  $annot_file1 = File::Spec->catfile($tmp, $file) if($file =~ m/^$index.$sample.*.1.ann$/);
	  $annot_file2 = File::Spec->catfile($tmp, $file) if($file =~ m/^$index.$sample.*.2.ann$/);
	}
	closedir($dh);

	my %gene_info;
	my $gene_gtf = $options->{'gene_gtf'};

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
	
	my %break1;
	open (my $ifh2, $annot_file1) or die "Could not open file '$annot_file1' $!";
	while (<$ifh2>) {
	  chomp;
		my $line1 = $_;
		my $break_annotation1 = parse_annotation($line1);
		my $gene1 = $break_annotation1->{'gene_name'};
		my $gene1_start = $gene_info{$gene1}{'feature_start'};
		my $gene1_end = $gene_info{$gene1}{'feature_end'};
		
		$break_annotation1->{'gene_start'} = $gene1_start;
		$break_annotation1->{'gene_end'} = $gene1_end;
		
		$break1{$break_annotation1->{'breakpoint'}} = $break_annotation1;
	}
	close ($ifh2);
		
	my %break2;
	open (my $ifh3, $annot_file2) or die "Could not open file '$annot_file2' $!";
	while (<$ifh3>) {
	  chomp;
		my $line2 = $_;
		
		my $break_annotation2 = parse_annotation($line2);
		my $gene2 = $break_annotation2->{'gene_name'};
		my $gene2_start = $gene_info{$gene2}{'feature_start'};
		my $gene2_end = $gene_info{$gene2}{'feature_end'};
		
		$break_annotation2->{'gene_start'} = $gene2_start;
		$break_annotation2->{'gene_end'} = $gene2_end;
		
		$break2{$break_annotation2->{'breakpoint'}} = $break_annotation2;
	}
	close ($ifh3);
	
	my $fusion;
	my $formatted_gene_line;
	my $formatted_exon_line;
	
	open(my $ofh1, '>', $gene_bedpe_file) or die "Could not open file $gene_bedpe_file $!";
	open(my $ofh2, '>', $exon_bedpe_file) or die "Could not open file $exon_bedpe_file $!";
	
	for my $brk (keys %break1){
	
	  $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
	    -breakpoint  => $brk,
	    -chr1	=> $break1{$brk}->{'chr'},
      -pos1_start	=> $break1{$brk}->{'pos_start'},
      -pos1_end	=> $break1{$brk}->{'pos_end'},
      -strand1	=> $break1{$brk}->{'strand'},
      -feature1	=> $break1{$brk}->{'feature'},
      -feature1_start	=> $break1{$brk}->{'feature_start'},
      -feature1_end	=> $break1{$brk}->{'feature_end'},
      -gene1	=> $break1{$brk}->{'gene_name'},
      -gene1_id	=> $break1{$brk}->{'gene_id'},
      -gene1_start	=> $break1{$brk}->{'gene_start'},
      -gene1_end	=> $break1{$brk}->{'gene_end'},
      -chr2	=> $break2{$brk}->{'chr'},
      -pos2_start	=> $break2{$brk}->{'pos_start'},
      -pos2_end	=> $break2{$brk}->{'pos_end'},
      -strand2	=> $break2{$brk}->{'strand'},
      -feature2	=> $break2{$brk}->{'feature'},
      -feature2_start	=> $break2{$brk}->{'feature_start'},
      -feature2_end	=> $break2{$brk}->{'feature_end'},
      -gene2	=> $break2{$brk}->{'gene_name'},
      -gene2_id	=> $break2{$brk}->{'gene_id'},
      -gene2_start	=> $break2{$brk}->{'gene_start'},
      -gene2_end	=> $break2{$brk}->{'gene_end'});

			if($break1{$brk}->{'feature'} eq "exon"){
				$fusion->exon1_num($break1{$brk}->{'exon_number'});
				$fusion->exon1_id($break1{$brk}->{'exon_id'});
				$fusion->transcript1_id($break1{$brk}->{'transcript_id'});
				$fusion->distance1($break1{$brk}->{'distance'});
			}

			if($break2{$brk}->{'feature'} eq "exon"){
				$fusion->exon2_num($break2{$brk}->{'exon_number'});
				$fusion->exon2_id($break2{$brk}->{'exon_id'});
				$fusion->transcript2_id($break2{$brk}->{'transcript_id'});
				$fusion->distance2($break2{$brk}->{'distance'});
			}
			
        $formatted_gene_line = $fusion->format_bedpe_line('gene');
        print $ofh1 $formatted_gene_line."\n";
        if(defined $fusion->exon1_num && $fusion->exon2_num){
          if($fusion->distance1 <= 10 && $fusion->distance2 <= 10){
            $formatted_exon_line = $fusion->format_bedpe_line('exon');
            print $ofh2 $formatted_exon_line."\n";
          }
        }
      
      
	  }
	close ($ofh1);
	close ($ofh2);
	
 	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
 	
  return 1;
}

sub filter_gtf {
  my ($gtf, $tmp, $feature) = @_;

  my $filtered_gtf = File::Spec->catfile($tmp, "filtered_$feature.gtf");
  
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
  
  return $filtered_gtf;
}

sub parse_annotation {
  my $line = shift;
  $line =~ s/"//g;
  
  my %annotation;
  my @fields = split "\t", $line;
  
  $annotation{'breakpoint'} = $fields[3];
  $annotation{'chr'} = $fields[0];
  $annotation{'pos_start'} = $fields[1];
  $annotation{'pos_end'} = $fields[2];
  $annotation{'strand'} = $fields[4];
  $annotation{'feature'} = $fields[7];
  $annotation{'feature_start'} = $fields[8];
  $annotation{'feature_end'} = $fields[9];
  
  my $annot_column = scalar @fields;
  $annotation{'distance'} = $fields[$annot_column-1];
	
	my @annot_fields = split /; /, $fields[$annot_column-2];
	foreach my $item(@annot_fields) {
  		my ($type,$value)= split / /, $item;
  		$annotation{$type} = $value;
	}

  return \%annotation;
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
  my ($line, $cols, $type) = @_;
  
  my @fields = split "\t", $line;
  
  my $row_length = scalar @fields;
  my $start = 0;
  $start = $row_length / 2  if($cols == 2);
  
  my $fusion = new Sanger::CGP::CompareFusions::FusionAnnotation(
	    -breakpoint  => $fields[$start + 7],
	    -chr1	=> $fields[$start],
      -strand1	=> $fields[$start + 8],
      -gene1	=> $fields[$start + 10],
      -gene1_id	=> $fields[$start + 11],
      -chr2	=> $fields[$start + 3],
      -strand2	=> $fields[$start + 9],
      -gene2	=> $fields[$start + 12],
      -gene2_id	=> $fields[$start + 13],
      -feature1	=> $fields[$start + 6]);
      
  if($type eq 'gene'){
		$fusion->gene1_start($fields[$start + 1]);
		$fusion->gene1_end($fields[$start + 2]);
		$fusion->gene2_start($fields[$start + 4]);
		$fusion->gene2_end($fields[$start + 5]);    
  }
      
  if($type eq 'exon'){
		$fusion->feature1_start($fields[$start + 1]);
		$fusion->feature1_end($fields[$start + 2]);
    $fusion->exon1_num($fields[$start + 14]);
		$fusion->exon1_id($fields[$start + 15]);
		$fusion->transcript1_id($fields[$start + 16]);
		$fusion->feature2_start($fields[$start + 4]);
		$fusion->feature2_end($fields[$start + 5]);
    $fusion->exon2_num($fields[$start + 17]);
    $fusion->exon2_id($fields[$start + 18]);
    $fusion->transcript2_id($fields[$start + 19]);
	}
  
  return $fusion;
}

sub process_exon_overlaps {
  my ($exon_overlap_file, $exon_list, $feature_source, $source1, $source2) = @_;
  
  # The fusion will be represented twice in the overlapping file, use the first set of columns by default
  my $col_set = 1;
  
  # If the first set of columns come from a defuse bedpe, use the second set of columns.
  # This is because star and tophat are consistent in terms of strand orientation. If we use
  # defuse data we would need to flip some of the gene and exon data
  $col_set = 2 if($source1 eq 'defuse');
  
  open (my $ifh, $exon_overlap_file) or die "Could not open file '$exon_overlap_file' $!";
	while (<$ifh>) {
	  chomp;
	  my $line = $_;
	  my @fields = split "\t", $line;
	  next if(($fields[15] ne $fields[37] && $fields[15] ne $fields[40]) || ($fields[18] ne $fields[37] && $fields[18] ne $fields[40]));
	  my $exon_fusion = parse_overlap($line, $col_set, 'exon');
	  $exon_list->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{$exon_fusion->{'breakpoint'}} = $exon_fusion if(!exists $exon_list->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{$exon_fusion->{'breakpoint'}});
	  my $brk1 = $fields[7];
	  my $brk2 = $fields[29];
	  $feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{'breakpoints'}{$brk1} = 1 if(!exists $feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{'breakpoints'}{$brk1});
	  $feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{'breakpoints'}{$brk2} = 1 if(!exists $feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{'breakpoints'}{$brk2});
		$feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{$source1} = 1 if(!exists $feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{$source1});
		$feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{$source2} = 1 if(!exists $feature_source->{$exon_fusion->{'exon1_id'}.':'.$exon_fusion->{'exon2_id'}}{$source2});
		
	}
	close ($ifh);
  
  return 1;
}

sub process_gene_overlaps {
  my ($gene_overlap_file, $gene_list, $feature_source, $source1, $source2) = @_;
  
  # The fusion will be represented twice in the overlapping file, use the first set of columns by default
  my $col_set = 1;
  
  # If the first set of columns come from a defuse bedpe, use the second set of columns.
  # This is because star and tophat are consistent in terms of strand orientation. If we use
  # defuse data we would need to flip some of the gene and exon data
  $col_set = 2 if($source1 eq 'defuse');
  
	open (my $ifh, $gene_overlap_file) or die "Could not open file '$gene_overlap_file' $!";
	while (<$ifh>) {
	  chomp;
	  my $line = $_;
	  my @fields = split "\t", $line;
	  next if(($fields[10] ne $fields[24] && $fields[10] ne $fields[26]) || ($fields[12] ne $fields[24] && $fields[12] ne $fields[26]));
		my $gene_fusion = parse_overlap($line, $col_set, 'gene');
	  $gene_list->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{$gene_fusion->{'breakpoint'}} = $gene_fusion if(!exists $gene_list->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{$gene_fusion->{'breakpoint'}});
	  my $brk1 = $fields[7];
	  my $brk2 = $fields[21];
	  $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{'breakpoints'}{$brk1} = 1 if(!exists $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{'breakpoints'}{$brk1});
		$feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{'breakpoints'}{$brk2} = 1 if(!exists $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{'breakpoints'}{$brk2});
	  $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{$source1} = 1 if(!exists $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{$source1});
	  $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{$source2} = 1 if(!exists $feature_source->{$gene_fusion->{'gene1_id'}.':'.$gene_fusion->{'gene2_id'}}{$source2});
	}
	close ($ifh);
  
  return 1;
}

sub run_bed_pairtopair {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

  my $sample = $options->{'sample'};
  
  my $prog = _which('bedtools');
  
  # There will always be at least two input files so build the command for the first comparison
  my $command1 = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, File::Spec->catfile($tmp, "1.$sample.gene.bedpe"),
  																										 File::Spec->catfile($tmp, "2.$sample.gene.bedpe"), 
  																										 File::Spec->catfile($tmp, "1_2.$sample.gene.bedpe_overlap");
  																										 
  my $command2 = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, File::Spec->catfile($tmp, "1.$sample.exon.bedpe"),
  																										 File::Spec->catfile($tmp, "2.$sample.exon.bedpe"), 
  																										 File::Spec->catfile($tmp, "1_2.$sample.exon.bedpe_overlap");
  my @commands = ($command1, $command2);
  
  if($options->{'num'} == 3){
    my $command3 = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, File::Spec->catfile($tmp, "1.$sample.gene.bedpe"),
  																										 	 File::Spec->catfile($tmp, "3.$sample.gene.bedpe"), 
  																										 	 File::Spec->catfile($tmp, "1_3.$sample.gene.bedpe_overlap");
  																										 
    my $command4 = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, File::Spec->catfile($tmp, "1.$sample.exon.bedpe"),
  																										 	 File::Spec->catfile($tmp, "3.$sample.exon.bedpe"),
  																											 File::Spec->catfile($tmp, "1_3.$sample.exon.bedpe_overlap");
  																											 
    my $command5 = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, File::Spec->catfile($tmp, "2.$sample.gene.bedpe"),
  																										 	 File::Spec->catfile($tmp, "3.$sample.gene.bedpe"),
  																										 	 File::Spec->catfile($tmp, "2_3.$sample.gene.bedpe_overlap");
  																										 
    my $command6 = $prog . sprintf $BEDTOOLS_PAIRTOPAIR, File::Spec->catfile($tmp, "2.$sample.exon.bedpe"),
  																										 	 File::Spec->catfile($tmp, "3.$sample.exon.bedpe"),
  																											 File::Spec->catfile($tmp, "2_3.$sample.exon.bedpe_overlap");																								 

	  push @commands, ($command3, $command4, $command5, $command6);
  }
  
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
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
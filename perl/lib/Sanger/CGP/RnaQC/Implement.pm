package Sanger::CGP::RnaQC::Implement;
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
use File::Which qw(which);
use FindBin qw($Bin);
use PCAP::Threaded;

use Sanger::CGP::CgpRna;

sub check_input {

  my $options = shift;
  
  my $inputdir = $options->{'indir'};
  my $sample = $options->{'sample'};
  
  # Check all the required input files exist
  PCAP::Cli::file_for_reading('bam-stats', File::Spec->catfile($inputdir, "$sample.bas"));
  PCAP::Cli::file_for_reading('bam-stats', File::Spec->catfile($inputdir, "$sample.transcriptome.bas"));
	PCAP::Cli::file_for_reading('gene-coverage', File::Spec->catfile($inputdir, "$sample.geneBodyCoverage.r"));
	PCAP::Cli::file_for_reading('rrna-percentage', File::Spec->catfile($inputdir, "$sample.rrna.txt"));
	PCAP::Cli::file_for_reading('read-distribution', File::Spec->catfile($inputdir, "$sample.read_dist.txt"));

  return 1;
}

sub gene_coverage {
  my $options = shift;
  
  my $inputdir = $options->{'indir'};
  my $sample = $options->{'sample'};
  my $gene_coverage_r_file = File::Spec->catfile($inputdir, "$sample.geneBodyCoverage.r");
  my $updated_coverage_r_file = File::Spec->catfile($inputdir, "$sample.geneBodyCoverage_UPDATED.r");
  
  # Modify the R script so that it runs the shapiro-wilk normality test and outputs that test statistic plus p-value to a text file. Also draw the gene body coverage graph with a more distinctly coloured line
  open (my $ifh, $gene_coverage_r_file) or die "Could not open file $gene_coverage_r_file $!";
  open (my $ofh, '>', $updated_coverage_r_file) or die "Could not open file $updated_coverage_r_file $!";
  
  while (<$ifh>){
    chomp;
    my $line = $_;
    if($line =~ m/^plot/){
      $line =~ s/col=icolor\[1\]/col="navy"/;
      print $ofh $line."\n";
    }
    else{
      print $ofh $line."\n";
    }
  }
  close($ifh);
  
  my $sample_r = $sample;
  my $output_file = File::Spec->catfile($inputdir, "$sample.gene.cov.bas");
  $sample_r =~ s/-/_/g;
  $sample_r = "V".$sample_r if($sample_r =~ m/^[0-9]/);
  print $ofh "test_result <- shapiro.test($sample_r)\n";
  print $ofh 'statistics<-unlist(test_result)'."\n";
  print $ofh 'shapiro_stat<-as.numeric(statistics["statistic.W"])'."\n";
  print $ofh 'shapiro_pval<-as.numeric(statistics["p.value"])'."\n";
  print $ofh 'output <- data.frame("gene_cov_stat"=shapiro_stat,"gene_cov_pval"=shapiro_pval)'."\n";
  print $ofh "write.table(output,file=\"$output_file\", sep=\"\t\", row.names=FALSE, col.names=TRUE, quote=FALSE)\n";

  close($ofh);
  
  my $command = _which('Rscript'); 
  $command .= " $updated_coverage_r_file";
  system($command);
  
  return 1;
}

sub junction_saturation {
  my $options = shift;
  my $inputdir = $options->{'indir'};
  my $sample = $options->{'sample'};
  my $saturation_r_file = File::Spec->catfile($inputdir, "$sample.junction_sat.junctionSaturation_plot.r");
  my $updated_saturation_r_file = File::Spec->catfile($inputdir, "$sample.junction_sat.junctionSaturation_plot_UPDATED.r");
  
  # Modify the R script so that it runs R commands to check the 'straightness' of the known junctions line. Capture a statistic for this into a bas file.
  open (my $ifh, $saturation_r_file) or die "Could not open file $saturation_r_file $!";
  open (my $ofh, '>', $updated_saturation_r_file) or die "Could not open file $updated_saturation_r_file $!";
  
  while (<$ifh>){
    chomp;
    my $line = $_;
    unless($line =~ m/^dev.off/){
      print $ofh $line."\n";
    }
  }
  close($ifh);
  
  my $output_file = File::Spec->catfile($inputdir, "$sample.junction.sat.bas");
  print $ofh "abline(mod <- lm(x ~ y))\n";
  print $ofh "gradient<-coef(mod)[2]\n";
  print $ofh 'output <- data.frame("junction_sat_stat"=gradient)'."\n";
  print $ofh "write.table(output,file=\"$output_file\", row.names=FALSE, col.names=TRUE, quote=FALSE)\n";
  close($ofh);
  
  my $command = _which('Rscript'); 
  $command .= " $updated_saturation_r_file";
  system($command);
  
  return 1;
}

sub read_distribution_stats {
  my $options = shift;
  
  my $inputdir = $options->{'indir'};
  my $sample = $options->{'sample'};
  
  my $read_dist_stats_file = File::Spec->catfile($inputdir, "$sample.read_dist.txt");
  my $read_dist_reformatted_file = File::Spec->catfile($inputdir, "$sample.read.dist.bas");
  my $stats_pattern = 'Exons|Introns|TSS_up_10kb|TES_down_10kb';
  my $total_reads = 0;
  my $exonic_reads = 0;
  my $intronic_reads = 0;
  my $intergenic_reads = 0;
  
  open (my $ifh, $read_dist_stats_file) or die "Could not open file $read_dist_stats_file $!";
  while (<$ifh>) {
		chomp;
		my $line = $_;
	  
	  if($line =~ m/Total Assigned Tags/){
	    my @fields = split " ", $line;
	    $total_reads = $fields[3];
	  }
	  
		if($line =~ m/$stats_pattern/){
			my @fields = split " ", $line;
			
			if($fields[0] =~ m/Exons/){
			  $exonic_reads += $fields[2];
			}
			elsif($fields[0] =~ m/Introns/){
			  $intronic_reads += $fields[2];
			}
			elsif($fields[0] =~ m/TSS|TES.*10kb/){
			  $intergenic_reads += $fields[2];
			}
		}
	}
  close($ifh);
  
  open(my $ofh, '>', $read_dist_reformatted_file) or die "Could not open file $read_dist_reformatted_file $!";
	
	print $ofh "#_total_read_dist_reads\t#_exonic_reads\t#_intronic_reads\t#_intergenic_within10kb_reads\n";
  print $ofh "$total_reads\t$exonic_reads\t$intronic_reads\t$intergenic_reads\n";
  
  close($ofh);
  
  return 1;
}

sub rrna_stats {
  my $options = shift;
  
  my $inputdir = $options->{'indir'};
  my $sample = $options->{'sample'};
  
  my $rrna_stats_file = File::Spec->catfile($inputdir, "$sample.rrna.txt");
  my $rrna_reformatted_file = File::Spec->catfile($inputdir, "$sample.rrna.bas");
  my $total_pattern = 'Total records';
  my $rrna_subset_pattern = "$sample.rRNA.in.bam";
  my $rrna_total_reads;
  my $rrna_subset_reads;
  
  open (my $ifh, $rrna_stats_file) or die "Could not open file $rrna_stats_file $!";
  while (<$ifh>) {
		chomp;
		my $line = $_;
		my @fields = split ":", $line;
		$fields[1] =~ s/\s+//;
		
		if($line =~ m/$total_pattern/){
			$rrna_total_reads = $fields[1];
		}
		elsif($line =~ m/$rrna_subset_pattern/){
		  $rrna_subset_reads = $fields[1];
		}
	}  
  close($ifh);
  
  open(my $ofh, '>', $rrna_reformatted_file) or die "Could not open file $rrna_reformatted_file $!";
  
  print $ofh "#_total_rrna_reads\t#_subset_rrna_reads\n";
  print $ofh "$rrna_total_reads\t$rrna_subset_reads\n";
  
  close($ofh);
  
  return 1;
}

sub transcriptome_stats {
  my $options = shift;
  
  my $inputdir = $options->{'indir'};
  my $sample = $options->{'sample'};
  
  my $first_bas_col_name = 'bam_filename';
  my $mean_insert_size_col_name = 'mean_insert_size';
  my $insert_size_sd_col_name = 'insert_size_sd';
  my $median_insert_size_col_name = 'median_insert_size';
  my $mean_insert_size_col_num;
  my $insert_size_sd_col_num;
  my $median_insert_size_col_num;
  my $mean_insert_size_val;
  my $insert_size_sd_val;
  my $median_insert_size_val;
  
  my $transcriptome_in_bas = File::Spec->catfile($inputdir, "$sample.transcriptome.bas");
  my $transcriptome_out_bas = File::Spec->catfile($inputdir, "$sample.insert.bas");
  
  open (my $ifh, $transcriptome_in_bas) or die "Could not open file $transcriptome_in_bas $!";
  while (<$ifh>) {
		chomp;
		my $line = $_;
		my @fields = split "\t", $line;
		my $length = scalar @fields;
		if($line =~ m/^$first_bas_col_name/){
		  for(my $i=0; $i < $length; $i++){
		    $mean_insert_size_col_num = $i if($fields[$i] eq $mean_insert_size_col_name);
		    $insert_size_sd_col_num = $i if($fields[$i] eq $insert_size_sd_col_name);
		    $median_insert_size_col_num = $i if($fields[$i] eq $median_insert_size_col_name);
		  }
		}
		else{
		  die "Transcriptome insert size values could not be found\n" if(!defined $mean_insert_size_col_num || !defined $insert_size_sd_col_num || !defined $median_insert_size_col_num);
		  $mean_insert_size_val = $fields[$mean_insert_size_col_num];
 			$insert_size_sd_val = $fields[$insert_size_sd_col_num];
  		$median_insert_size_val = $fields[$median_insert_size_col_num];
		}
	}
  close($ifh);
  
  open(my $ofh, '>', $transcriptome_out_bas) or die "Could not open file $transcriptome_out_bas $!";
  print $ofh "transcriptome_$mean_insert_size_col_name\ttranscriptome_$insert_size_sd_col_name\ttranscriptome_$median_insert_size_col_name\n";
  print $ofh "$mean_insert_size_val\t$insert_size_sd_val\t$median_insert_size_val\n";
  close($ofh);
  
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
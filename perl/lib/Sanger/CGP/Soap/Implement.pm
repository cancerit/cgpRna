package Sanger::CGP::Soap::Implement;
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
use File::Spec;
use File::Path qw(remove_tree make_path);
use File::Basename;
use File::Copy;
use Cwd;
use FindBin qw($Bin);
use Capture::Tiny qw(capture capture_stdout);
use PCAP::Cli;
use PCAP::Threaded;
use PCAP::Bwa::Meta;
use PCAP::Bam;
use Sanger::CGP::Soap;

const my $BAMFASTQ => q{ exclude=QCFAIL,SECONDARY,SUPPLEMENTARY T=%s S=%s O=%s O2=%s gz=1 level=1 F=%s F2=%s filename=%s};
const my $FUSIONS_FILTER => q{ -i %s -s %s -n %s -o %s -p soap};
const my $SOAP_MAX_CORES => 8;
const my $SOAP => q{ %s -c %s -fm %d -o %s -tp %s -l %s -fd %s };

sub check_input {
	my $options = shift;

	my $ref_data = $options->{'refdataloc'};
	my $species = $options->{'species'};
	my $ref_build = $options->{'referencebuild'};
	my $gene_build = $options->{'genebuild'};
	my $ref_build_loc = File::Spec->catdir($ref_data, $species, $ref_build);

	# Check the normal fusions file exists
	PCAP::Cli::file_for_reading('normals-list',File::Spec->catfile($ref_build_loc,$options->{'normalfusionslist'}));
	
	my $input_meta = PCAP::Bwa::Meta::files_to_meta($options->{'tmp'}, $options->{'raw_files'}, $options->{'sample'});
	
	$options->{'meta_set'} = $input_meta;
	
	# Check the input data.
	$options->{'max_split'} = scalar @{$input_meta};
	
	my $fq = 0;
	my $fqgz = 0;
	my $bam = 0;
	
	# If the input is fastq need to check whether there is a mixture of gzipped and uncompressed fastqs.
	
	for my $input(@{$input_meta}) {
		if($input->fastq) {
			if($input->fastq =~ m/\.gz$/) {
				$fqgz = 1;
			}
			else{
				$fq = 1;
			}
		}
		else{
			$bam = 1;
		}
	}
	
	# If BAMs have been detected set flag so that they will be converted to pairs of fastq files in the prepare subroutine
	if($bam == 1){
		$options->{'bam'} = 1;
	}
	else{
		# If a mixture of gzipped and uncompressed fastqs, set flag to gzip the plain fastq/fq files in the prepare subroutine
		if($fq == $fqgz){
			$options->{'mixedfq'} = 1;
		}
	}
	
	return 1;
}

sub filter_fusions {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $sample = $options->{'sample'};
	my $soap_outdir = File::Spec->catdir($options->{'tmp'}, 'soapfuse', 'final_fusion_genes', $sample);
	my $fusions_file = File::Spec->catfile($soap_outdir, "$sample.final.Fusion.specific.for.genes");
	die "Please run the soap step prior to filter\n" unless(-d $soap_outdir);
	die "One of the SOAPfuse output files is missing, please run the soap step prior to filter.\n" unless(-e $fusions_file);

	my $normals_file = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},$options->{'normalfusionslist'});

	my $command = "$^X ";
	$command .= _which('filter_fusions.pl');
	$command .= sprintf $FUSIONS_FILTER, 	$fusions_file,
						$sample,
						$normals_file,
						$options->{'outdir'};

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	return 1;
}

sub prepare {

	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $sample = $options->{'sample'};
 	my $library = $options->{'library'};
	my $inputdir = File::Spec->catdir($tmp, 'input', $sample, $library);
	
	if($options->{'bam'} || $options->{'mixedfq'}){
		my $input_meta = $options->{'meta_set'};
		my $iter = 1;
		
		for my $input(@{$input_meta}) {
			next if($iter++ != $index); # skip to the relevant element in the list
			# If the file name ends fastq or fq we need to gzip it to make it consistent with other input files in the list
			if($input->fastq) {
				if($input->fastq =~ m/f(ast)?q$/) {
					my ($filename, $filepath) = fileparse($input->in);
					my $suffix = ".".$input->fastq;
					my $outfile = File::Spec->catfile($inputdir,$filename);
					my $raw_files = $options->{'raw_files'};
					
					if($input->paired_fq){

						my $infile1 = $input->in."_1".$suffix;
						my $infile2 = $input->in."_2".$suffix;
						my $outfile1 = $outfile."_1".$suffix.".gz";
						my $outfile2 = $outfile."_2".$suffix.".gz";
						
						system([0,2], "(gzip -c $infile1 > $outfile1) >& /dev/null") unless(-e $outfile1);
						system([0,2], "(gzip -c $infile2 > $outfile2) >& /dev/null") unless(-e $outfile2);
						
						for(my $i=0; $i<@$raw_files; $i++){
							@$raw_files[$i] = $outfile1 if(@$raw_files[$i] =~ $infile1);
							@$raw_files[$i] = $outfile2 if(@$raw_files[$i] =~ $infile2);
						}
					}
					else{

						my $infile = $input->in.$suffix;
						my $outfile = $outfile.$suffix.".gz";
						
						system([0,2], "(gzip -c $infile > $outfile) >& /dev/null") unless(-e $outfile);
						
						for(my $i=0; $i<@$raw_files; $i++){
							@$raw_files[$i] = $outfile if(@$raw_files[$i] =~ $infile);
						}				
					}
				}
			}
			else{
				# It's a BAM so call bamtofastq to split into gzipped paired fastq files
				my $rg = $input->{'rg'};
				
				my $command = _which('bamtofastq') || die "Unable to find 'bamtofastq' in path";
				$command .= sprintf $BAMFASTQ, 	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.s"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o1"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o2"),
								File::Spec->catfile($inputdir, $sample.'.'.$rg.'_1.fastq.gz'),
								File::Spec->catfile($inputdir, $sample.'.'.$rg.'_2.fastq.gz'),
								$input->in;
																	
					PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
			
			}
		}
	}
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
	
	return 1;
}

sub soap {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $threads = $SOAP_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $SOAP_MAX_CORES);
	my $sample = $options->{'sample'};
	my $soap = $options->{'soappath'};
	my $inputdir = File::Spec->catdir($tmp, 'input');
	my $outdir = File::Spec->catdir($tmp, 'soapfuse');
	my $sample_file = File::Spec->catfile($tmp, 'sample-list');
	die "$sample_file could not be found in the temporary directory, please check and run soapsetup before the soap stage\n" unless(-e $sample_file);
	
	# Get the relevant defuse config file for the reference and gene builds
	my $soap_config = File::Spec->catfile($options->{'refdataloc'}, $options->{'species'}, $options->{'referencebuild'}, $options->{'genebuild'}, $options->{'soapconfig'} );
	my $command = sprintf $SOAP,	$soap,
					$soap_config,
					$threads,
					$outdir,
					$sample,
					$sample_file,
					$inputdir;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  
  return 1;
}

sub soap_setup {
  my $options = shift;
  
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);	  

  # Check the input files are in the right location and write the sample list file to the temporary directory
  my $sample = $options->{'sample'};
 	my $library = $options->{'library'};
 	my $readlength = $options->{'readlength'};
	my $inputdir = File::Spec->catdir($tmp, 'input', $sample, $library);
	my $sample_list_file = File::Spec->catfile($tmp, 'sample-list');
	
	my @input1;
  my @input2;
  open(my $ofh1, '>', $sample_list_file) or die "Could not open file $sample_list_file $!";
  my $input_meta = $options->{'meta_set'};
  for my $input(@{$input_meta}) {
    if($input->fastq) {
			# Paired fastq input
			if($input->paired_fq) {
			  # TODO
			}
			# Interleaved fastq input
			else{
				# TODO
			}
		}
  	else{
			my $rg = $input->rg;
			die "Please run the prepare step prior to soapsetup\n" unless(-e File::Spec->catfile($inputdir, $sample.'.'.$rg.'_1.fastq.gz'));
			my $filename = $sample.'.'.$rg;
		  print $ofh1 "$sample\t$library\t$filename\t$readlength\n";
		}
  }
  close($ofh1);
  
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

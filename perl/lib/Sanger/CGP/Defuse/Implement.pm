package Sanger::CGP::Defuse::Implement;
##########LICENCE ##########
#Copyright (c) 2015 Genome Research Ltd.
###
#Author: Cancer Genome Project <cgpit@sanger.ac.uk>
###
#This file is part of cgpRna.
###
#TopHatFusion is free software: you can redistribute it and/or modify it under
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
use Sanger::CGP::Defuse;

use Data::Dumper;

const my $DEFUSE => q{ %s -c %s -p %d -o %s -n %s -l %s -s direct -1 %s -2 %s};
const my $BAMFASTQ => q{ exclude=QCFAIL,SECONDARY,SUPPLEMENTARY T=%s S=%s O=%s O2=%s level=1 F=%s F2=%s filename=%s};
const my $FUSIONS_FILTER => q{ -i %s -s %s -n %s -o %s};
const my $DEFUSE_MAX_CORES => 16;

sub check_input {
	my $options = shift;
	
	my $ref_data = $options->{'refdataloc'};
	my $species = $options->{'species'};
	my $ref_build = $options->{'referencebuild'};
	my $gene_build = $options->{'genebuild'};
	my $ref_build_loc = File::Spec->catdir($ref_data, $species, $ref_build);

	# Check the normal fusions file exists for the filtering step
	PCAP::Cli::file_for_reading('normals-list',File::Spec->catfile($ref_build_loc,'normal-fusions',$options->{'normalfusionslist'}));
	
	# Check the defuse config file exists for the ref-gene build
	PCAP::Cli::file_for_reading('defuse-config',File::Spec->catfile($ref_build_loc,$gene_build,$options->{'defuseconfig'}));
	
	$options->{'meta_set'} = PCAP::Bwa::Meta::files_to_meta($options->{'tmp'}, $options->{'raw_files'}, $options->{'sample'});	
	
	# Check the input data type
	my $input_meta = $options->{'meta_set'};
	my $input = @{$input_meta}[0];
	
	# If an interleaved fastq has been entered, stop the script as this is not valid input for deFuse which only takes paired fastq files
	die "Interleaved or single-end fastqs are not valid for deFuse, please re-try with BAM or paired fastq input" if($input->fastq && !$input->paired_fq);
		
	# If the input includes BAM files update flag in options to trigger the bamtofastq subroutine
	$options->{'bam'} = 1 unless($input->fastq);
	$options->{'max_split'} = scalar @{$options->{'meta_set'}};
		
	return 1;
}

sub compress_sam {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my $sample = $options->{'sample'};
	my $defuse_outdir = File::Spec->catdir($options->{'tmp'}, "defuse_$sample");
	my $in_sam = File::Spec->catfile($defuse_outdir, 'cdna.pair.sam');
	my $sam_gz = File::Spec->catfile($options->{'outdir'}, 'cdna.pair.sam.gz');	
	
	my $command = _which('gzip');
	$command .= sprintf ' -c %s > %s', $in_sam, $sam_gz;
	
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);	
	return 1;
}

sub defuse {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $threads = $DEFUSE_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $DEFUSE_MAX_CORES);
	my $sample = $options->{'sample'};
	my $defuse = $options->{'defusepath'};
	my $outdir = File::Spec->catdir($tmp, "defuse_$sample");
	my $fastq1;
	my $fastq2;
	
	# If the input was BAM or multiple files then the paired fastqs will reside in the /tmp/input folder, otherwise grab the pair from the raw_files array on the options hash
	if($options->{'bam'} || $options->{'max_split'} > 1){
		$fastq1 = File::Spec->catfile($tmp,'input', $sample."_1.fastq");
		$fastq2 = File::Spec->catfile($tmp,'input', $sample."_2.fastq");
	}
	else {
		my $raw_files = $options->{'raw_files'};
		my @sorted_files = sort @{$raw_files};
		$fastq1 = $sorted_files[0];
		$fastq2 = $sorted_files[1];
	}
	
	# Get the relevant defuse config file for the reference and gene builds
	my $defuse_config = File::Spec->catfile($options->{'refdataloc'}, $options->{'species'}, $options->{'referencebuild'}, $options->{'genebuild'}, $options->{'defuseconfig'} );
	my $command = sprintf $DEFUSE,	$defuse,
																	$defuse_config,
																	$threads,
																	$outdir,
																	$sample,
																	$outdir,
																	$fastq1,
																	$fastq2;
																	
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
	return 1;
}

sub filter_fusions {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $sample = $options->{'sample'};
	my $defuse_outdir = File::Spec->catdir($options->{'tmp'}, "defuse_$sample");
	my $fusions_file = File::Spec->catfile($defuse_outdir, 'results.filtered.tsv');
	die "Please run the defuse step prior to filter\n" unless(-d $defuse_outdir);
	die "One of the deFuse output files is missing, please run the defuse step prior to filter.\n" unless(-e $fusions_file && -e File::Spec->catfile($defuse_outdir, 'cdna.pair.sam'));

	my $normals_file = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},'normal-fusions',$options->{'normalfusionslist'});

	my $command = "$^X ";
	$command .= _which('filter_defuse_fusions.pl');
	$command .= sprintf $FUSIONS_FILTER, 	$fusions_file,
																				$sample,
																				$normals_file,
																				$options->{'outdir'};

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	return 1;
}

sub merge {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $inputdir = File::Spec->catdir($tmp, 'input');
	my $sample = $options->{'sample'};
	my @files1;
	my @files2;
	my @commands;
	
	# If the input is BAM then all input fastqs will reside in the ../tmp/input directory so can search for _1 and _2 to find the files to merge
	if($options->{'bam'}){
 		opendir(my $dh, $inputdir);
		while(my $file = readdir $dh) {
			push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/_1.f/);
			push @files2, File::Spec->catfile($inputdir, $file) if($file =~ m/_2.f/);
		}
		closedir($dh);
		
		my $infiles1 = join(' ', sort @files1);
		my $infiles2 = join(' ', sort @files2);
		my $outfile1 = File::Spec->catfile($inputdir,$sample."_1.fastq");
		my $outfile2 = File::Spec->catfile($inputdir,$sample."_2.fastq");
		
		push @commands, "cat $infiles1 > $outfile1";
		push @commands, "cat $infiles2 > $outfile2";

	}
	else {

		my $raw_files = $options->{'raw_files'};
		my $input_meta = $options->{'meta_set'};
		my $input = @{$input_meta}[0];
		
		for my $file(@{$raw_files}) {
			push @files1, $file if($file =~ m/_1.f/);
			push @files2, $file if($file =~ m/_2.f/);
		}
			
		my $infiles1 = join(' ', sort @files1);
		my $infiles2 = join(' ', sort @files2);
		my $outfile1 = File::Spec->catfile($inputdir,$sample."_1.fastq");
		my $outfile2 = File::Spec->catfile($inputdir,$sample."_2.fastq");
			
		push @commands, "cat $infiles1 > $outfile1";
		push @commands, "cat $infiles2 > $outfile2";
	}
	
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

	return 1;
}

sub prepare {

	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	my $input_meta = $options->{'meta_set'};
	my $iter = 1;
	for my $input(@{$input_meta}) {
		next if($iter++ != $index); # skip to the relevant element in the list
		# If the file is a gzipped fastq or fq we need to decompress it for deFuse to accept it
		if($input->fastq) {
			if($input->fastq =~ m/\.gz$/) {
				my $inputdir = File::Spec->catdir($tmp, 'input');
				my ($filename, $filepath) = fileparse($input->in);
				my $suffix = $input->fastq;
				$suffix =~ s/.gz//;
				my $outfile = File::Spec->catfile($inputdir,$filename);
					 
				my $infile1 = $input->in."_1.".$input->fastq;
				my $infile2 = $input->in."_2.".$input->fastq;
				my $outfile1 = $outfile."_1";
				my $outfile2 = $outfile."_2";
					 
				system([0,2], "(gunzip -c $infile1 > $outfile1.$suffix) >& /dev/null") unless(-e "$outfile1.$suffix");
				system([0,2], "(gunzip -c $infile2 > $outfile2.$suffix) >& /dev/null") unless(-e "$outfile2.$suffix");
						
				my $raw_files = $options->{'raw_files'};
				for(my $i=0; $i<@$raw_files; $i++){
					@$raw_files[$i] = $outfile1.".".$suffix if(@$raw_files[$i] =~ $infile1);
					@$raw_files[$i] = $outfile2.".".$suffix if(@$raw_files[$i] =~ $infile2);
				}
			}
		}
		else{
			my $inputdir = File::Spec->catdir($tmp, 'input');
			my $sample = $options->{'sample'};
			my $rg = $input->{'rg'};
			my $command = _which('bamtofastq') || die "Unable to find 'bamtofastq' in path";
			$command .= sprintf $BAMFASTQ, File::Spec->catfile($tmp, "bamtofastq.$sample.$rg"),
																	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.s"),
																	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o1"),
																	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o2"),
																	File::Spec->catfile($inputdir, $sample.'.'.$rg.'_1.fastq'),
																	File::Spec->catfile($inputdir, $sample.'.'.$rg.'_2.fastq'),
																	$input->in;
																	
			PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
		}
	}
	
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
	
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

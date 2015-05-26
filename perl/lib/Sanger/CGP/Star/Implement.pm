package Sanger::CGP::Star::Implement;
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
use Sanger::CGP::Star;

use Data::Dumper;

const my $BAMFASTQ => q{ exclude=QCFAIL,SECONDARY,SUPPLEMENTARY T=%s S=%s O=%s O2=%s gz=1 level=1 F=%s F2=%s filename=%s};
const my $FUSIONS_FILTER => q{ -i %s -s %s -n %s -o %s -p star};
const my $STAR_MAX_CORES => 16;
const my $STAR_DEFAULTS_SECTION => 'star-parameters';
const my $STAR_FUSION_SECTION => 'star-fusion-parameters';
const my $STAR => q{ %s %s --readFilesIn %s };
const my $STAR_FUSION => q{ %s --chimeric_out_sam %s --chimeric_junction %s --ref_GTF %s --min_novel_junction_support 10 --min_alt_pct_junction 10.0 --out_prefix %s };
const my $SAMTOBAM => q{ view -bS %s > %s };


sub check_input {
	my $options = shift;

	my $ref_data = $options->{'refdataloc'};
	my $species = $options->{'species'};
	my $ref_build = $options->{'referencebuild'};
	my $gene_build = $options->{'genebuild'};
	my $ref_build_loc = File::Spec->catdir($ref_data, $species, $ref_build);

	# Check the gtf and normal fusions files exist
	PCAP::Cli::file_for_reading('gtf-file', File::Spec->catfile($ref_build_loc, $gene_build, $options->{'gtffilename'}));
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
	my $star_outdir = File::Spec->catdir($options->{'tmp'}, 'star');
	my $fusions_file = File::Spec->catfile($star_outdir, "$sample.fusion_candidates.txt");
	die "The star fusion output files are missing, please run the starfusion step prior to filter.\n" unless(-e $fusions_file);

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
	
	if($options->{'bam'} || $options->{'mixedfq'}){
		my $input_meta = $options->{'meta_set'};
		my $inputdir = File::Spec->catdir($tmp, 'input');
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
				my $sample = $options->{'sample'};
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

sub process_star_params {
	my ($options, $fusion_mode) = @_;
	
	my $ini_file = $options->{'config'};
	my $cfg = new Config::IniFiles( -file => $ini_file ) or die "Could not open config file: $ini_file";
	
	die "The star default parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($STAR_DEFAULTS_SECTION));
	
	# Get the Star sample specific parameters from the options hash
	my $threads = $STAR_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $STAR_MAX_CORES);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'runThreadN', $threads);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'outFileNamePrefix', $options->{'tmp'}.'/star/');
	
	my $gtf = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'}, $options->{'genebuild'}, $options->{'gtffilename'});
	my $star_index = File::Spec->catdir($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'}, 'star-index');
	
	$cfg->setval($STAR_DEFAULTS_SECTION, 'sjdbGTFfile', $gtf);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'genomeDir', $star_index);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'outSAMattrRGline', $options->{'rgline'}) if(defined $options->{'rgline'});
	
	my @star_command;
	my @star_defaults = $cfg->Parameters($STAR_DEFAULTS_SECTION);
	
	for my $key(@star_defaults){
		if($cfg->val($STAR_DEFAULTS_SECTION, $key) ne ''){
			push @star_command, "--".$key." ".$cfg->val($STAR_DEFAULTS_SECTION, $key);
		}
	}
	
	if(defined $fusion_mode) {
		die "The star-fusion default parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($STAR_FUSION_SECTION));
		
		my @star_fusion_defaults = $cfg->Parameters($STAR_FUSION_SECTION);
		for my $key(@star_fusion_defaults){
			if($cfg->val($STAR_FUSION_SECTION, $key) ne ''){
				push @star_command, "--".$key." ".$cfg->val($STAR_FUSION_SECTION, $key);
			}
		}
		
	}
	
	return join(" ", @star_command);
}

sub prog_version {
	my $options = shift;
	my $star_path = $options->{'starpath'};

	if(! defined $options->{'starpath'} || $options->{'starpath'} eq ''){
		$star_path = _which('STAR');
		$options->{'starath'} = $star_path;
	}

	my $star_version;
	{
		no autodie qw(system);
		my ($stdout, $stderr, $exit) = capture{ system("$star_path --version"); };
		($star_version) = $stdout =~ /STAR_([[:digit:]\.]+.*)/m;
	}
	return $star_version;
}

sub sam_to_bam {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $sample = $options->{'sample'};
	my $outdir = $options->{'outdir'};
	my $star_outdir = File::Spec->catdir($options->{'tmp'}, 'star');
	my $command .= _which('samtools');
	$command .= sprintf $SAMTOBAM,	File::Spec->catfile($star_outdir, 'Chimeric.out.sam'),
																	File::Spec->catfile($outdir, $sample.'.Chimeric.out.bam');

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
	return 1;
}

sub star_chimeric {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $threads = $STAR_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $STAR_MAX_CORES);
	my $sample = $options->{'sample'};
	
	# Check the input data
	my $input_meta = $options->{'meta_set'};
	
	# Get the RG header information to format parameters --rg-id and --rg-sample
	my $first_file = $input_meta->[0];
	my $rg_line;

	if($first_file->fastq) {
		$rg_line = $first_file->rg_header(q{\t});
	}
	else {
		($rg_line, undef) = PCAP::Bam::rg_line_for_output($first_file->in, $sample);
		$rg_line = $rg_line;
	}
	
	# Format the @RG header line for the output BAM file. Quotes need to be around the description (DS:) tag text
	$rg_line =~ s/^(.*)(DS:[^\\]+)(\\t.*$)/$1"$2"$3/;
	$rg_line =~ s/^\@RG\\t//;
	$rg_line =~ s/\\t/ /g;
	$options->{'rgline'} = $rg_line;
	
	my $star_params = process_star_params($options, 1);
	
	my @files1;
	my @files2;
	my $infiles1;
	my $infiles2;

	# If the input is BAM then all input fastqs will reside in the ../tmp/input directory so can search for _1 and _2 to find the files
	if($options->{'bam'}){
		my $inputdir = File::Spec->catdir($tmp, 'input');
 		opendir(my $dh, $inputdir);
		while(my $file = readdir $dh) {
			push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/_1.f/);
			push @files2, File::Spec->catfile($inputdir, $file) if($file =~ m/_2.f/);
		}
		closedir($dh);
		
		$infiles1 = join(',', sort @files1);
		$infiles2 = join(',', sort @files2);
	}
	else{
		my $raw_files = $options->{'raw_files'};
		
		# If there are multiple input fastqs, some gzipped and some not, need to check both the input directory and raw_files array for the locations of the input files
		if($options->{'mixedfq'}){
			my $inputdir = File::Spec->catdir($tmp, 'input');
			if($first_file->paired_fq){
				for my $file(@{$raw_files}) {
					push @files1, $file if($file =~ m/_1.f.*q.gz$/);
					push @files2, $file if($file =~ m/_2.f.*q.gz$/);
				}
				opendir(my $dh, $inputdir);
				while(my $file = readdir $dh) {
					push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/_1.f.*q.gz$/);
					push @files2, File::Spec->catfile($inputdir, $file) if($file =~ m/_2.f.*q.gz$/);
				}
				closedir($dh);
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = join(',', sort @files2);
			
			}
			else{
				for my $file(@{$raw_files}) {
					push @files1, $file if($file =~ m/.f.*q.gz$/);
				}
				opendir(my $dh, $inputdir);
				while(my $file = readdir $dh) {
					push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/.f.*q.gz$/);
				}
				closedir($dh);
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = '';
			}

		}
		# Otherwise there should be just one pair of fastq files or an interleaved/single ended fastq file so only need to check the raw_files array
		else{
			if($first_file->paired_fq){
				for my $file(@{$raw_files}) {
					push @files1, $file if($file =~ m/_1.f/);
					push @files2, $file if($file =~ m/_2.f/);
				}
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = join(',', sort @files2);
				
			}
			else{
				for my $file(@{$raw_files}) {
					push @files1, $file;
				}
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = '';
				
			}
		}
	}
	
	my $command = sprintf $STAR,	$options->{'starpath'},
					$star_params,
					$infiles1." ".$infiles2;
																
	$command .= "--readFilesCommand zcat" if($infiles1 =~ m/\.gz$/);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
	return 1;
}

sub star_fusion {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $star_dir = File::Spec->catdir($tmp, 'star');
	my $chimeric_junction = File::Spec->catfile($star_dir, 'Chimeric.out.junction');
	my $chimeric_sam = File::Spec->catfile($star_dir, "Chimeric.out.sam");
	die "Some star-fusion setup files are missing, please check the output from STAR and re-run the previous stage (star) if necessary." unless(-e $chimeric_junction && -e $chimeric_sam);
	
	my $sample = $options->{'sample'};
	my $gtf = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'}, $options->{'genebuild'}, $options->{'gtffilename'});
	
	my $command = sprintf $STAR_FUSION,	$options->{'starfusionpath'},
						$chimeric_sam,
						$chimeric_junction,
						$gtf,
						$star_dir."/$sample";
																			
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
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

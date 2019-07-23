package Sanger::CGP::Tophat::Implement;
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
use Sanger::CGP::CgpRna;

const my @BOWTIE1_SUFFIXES => qw(.1.ebwt .2.ebwt .3.ebwt .4.ebwt .rev.1.ebwt .rev.2.ebwt .fa .fa.fai);
const my @BOWTIE2_SUFFIXES => qw(.1.bt2 .2.bt2 .3.bt2 .4.bt2 .rev.1.bt2 .rev.2.bt2 .fa .fa.fai);
const my $BAMFASTQ => q{ exclude=SECONDARY,SUPPLEMENTARY T=%s S=%s O=%s O2=%s gz=1 level=1 F=%s F2=%s filename=%s};
const my $FUSIONS_FILTER => q{ -i %s -s %s -n %s -o %s -p tophat};
const my $ADD_STRAND => q{ -i %s -s %s -p %s -o %s};
const my $FUSIONS_SPLIT => 10000;
const my $TOPHAT_MAX_CORES => 16;
const my $TOPHAT_DEFAULTS_SECTION => 'tophat-parameters';
const my $TOPHAT_FUSION_SECTION => 'tophat-fusion-parameters';
const my $TOPHAT_POST_SECTION => 'tophat-post-parameters';

sub add_strand {
	my $options = shift;

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my $sample = $options->{'sample'};
	my $post_outdir = File::Spec->catdir($options->{'tmp'}, 'tophatpostrun/tophatfusion_'.$sample);
	my $fusions_file = File::Spec->catfile($post_outdir, "$sample.tophat-fusion.normals.filtered.txt");
	my $potential_fusions = File::Spec->catfile($post_outdir, 'potential_fusion.txt');
	die "Please run the filter_fusions step prior to filter\n" unless(-e $fusions_file && -e $potential_fusions);

	if(-s $potential_fusions == 0) {
		system("echo '##EOF##' > $potential_fusions") && die "An error occurred: $!";
	}

	my $command = "$^X ";
	$command .= _which('tophat_add_strand.pl');
	$command .= sprintf $ADD_STRAND, 	$fusions_file,
						$sample,
						$potential_fusions,
						$post_outdir;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

	return 1;
}

sub bam_to_fastq {

	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);

	my $command;
	my $rg;
	my $sample = $options->{'sample'};
	my $inputdir = File::Spec->catdir($tmp, 'input');
	my $input_meta = $options->{'meta_set'};
	my $iter = 1;
	for my $input(@{$input_meta}) {
		next if($iter++ != $index); # skip to the relevant element in the list
		$command = _which('bamtofastq') || die "Unable to find 'bamtofastq' in path";
		$rg = $input->rg;
		$command .= sprintf $BAMFASTQ, 	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg"),
						File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.s"),
						File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o1"),
						File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o2"),
						File::Spec->catfile($inputdir, $sample.'.'.$rg.'_1.fastq.gz'),
						File::Spec->catfile($inputdir, $sample.'.'.$rg.'_2.fastq.gz'),
						$input->in;
	}

 	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);

	return 1;
}

sub check_input {
	my $options = shift;

	# Check the bowtie reference index files are in place and are appropriate for the version of bowtie being used
	my $suffix;
	my $ref_prefix = $options->{'referenceindex'};
	my $trans_prefix = $options->{'transcriptomeindex'};
	my $refdata = File::Spec->catdir($options->{'refdataloc'},$options->{'species'});
	my $ens_refdata = File::Spec->catdir($refdata,$options->{'referencebuild'});

	if($options->{'bowtieversion'} == 1){
		for $suffix(@BOWTIE1_SUFFIXES){
			PCAP::Cli::file_for_reading('bowtie1-ref-index',File::Spec->catfile($ens_refdata,$ref_prefix.$suffix));
			PCAP::Cli::file_for_reading('bowtie1-transcriptome-index',File::Spec->catfile($ens_refdata,'tophat',$options->{'genebuild'},$trans_prefix.$suffix));
			$options->{'referencepath'} = File::Spec->catfile($ens_refdata,$ref_prefix);
			$options->{'transcriptomepath'} = File::Spec->catfile($ens_refdata,'tophat',$options->{'genebuild'},$trans_prefix);
		}
	}
	else{
		for $suffix(@BOWTIE2_SUFFIXES){
			PCAP::Cli::file_for_reading('bowtie2-ref-index',File::Spec->catfile($ens_refdata,$ref_prefix.$suffix));
			PCAP::Cli::file_for_reading('bowtie2-transcriptome-index',File::Spec->catfile($ens_refdata,'tophat',$options->{'genebuild'},$trans_prefix.$suffix));
			$options->{'referencepath'} = File::Spec->catfile($ens_refdata,$ref_prefix);
			$options->{'transcriptomepath'} = File::Spec->catfile($ens_refdata,'tophat',$options->{'genebuild'},$trans_prefix);
		}
	}

	# Check the TopHat Fusion Post files exist
	my $ucsc_prefix = $options->{'tophatpostindex'};
	for $suffix(@BOWTIE1_SUFFIXES){
		PCAP::Cli::file_for_reading('bowtie1-tophatpost-index',File::Spec->catfile($ens_refdata,'tophat',$ucsc_prefix.$suffix));
	}
	$options->{'tophatpostpath'} = File::Spec->catfile($ens_refdata,'tophat',$ucsc_prefix);
	PCAP::Cli::file_for_reading('refGene',File::Spec->catfile($ens_refdata,'tophat',$options->{'refgene'}));
	PCAP::Cli::file_for_reading('ensGene',File::Spec->catfile($ens_refdata,'tophat',$options->{'genebuild'},$options->{'ensgene'}));

	# Check the normal fusions file exists for the filtering step
	PCAP::Cli::file_for_reading('normals-list',File::Spec->catfile($ens_refdata,'cgpRna',$options->{'normalfusionslist'}));

	$options->{'meta_set'} = PCAP::Bwa::Meta::files_to_meta($options->{'tmp'}, $options->{'raw_files'}, $options->{'sample'});

	# If the input includes BAM files update flag in options to trigger the bamtofastq subroutine
	my $input_meta = $options->{'meta_set'};
	# Only need to check the first input file to see if it's BAM
	my $input = @{$input_meta}[0];
	$options->{'bam'} = 1 unless($input->fastq);
	$options->{'max_split'} = scalar @{$options->{'meta_set'}};

	return 1;
}

sub check_ref_seqs {
	my $options = shift;

	my $tophat_ref = $options->{'referencepath'}.'.fa.fai';
	my $tophatpost_ref = $options->{'tophatpostpath'}.'.fa.fai';

	my $tophat_fai_seqs = capture_stdout { system('cut', '-f', 1, $tophat_ref ); };
	my $tophatpost_fai_seqs = capture_stdout { system('cut', '-f', 1, $tophatpost_ref ); };

	my @tophat_all_seqs = split /\n/, $tophat_fai_seqs;
	my @tophatpost_all_seqs = split /\n/, $tophatpost_fai_seqs;

	# The current assumption is that tophatpost always needs the chromosome names to be prefixed with chr.
	# When comparing, keep the reference chromosome names as is and update the tophatpost names accordingly
	# that way we have the correctly formatted exclusion string ready for the tophat fusion mapping parameter
	if($tophat_all_seqs[0] !~ /^chr/){
		for (@tophatpost_all_seqs){
			s/^chr//;
		}
	}

	# First cross compare the tophat reference with the tophat post UCSC reference to see what's common
	my %ref = ();
	my @isect = ();

	map { $ref{$_} = 1 } @tophat_all_seqs;
	@isect = grep { $ref{$_} } @tophatpost_all_seqs;

	# Subtract the common chromosomes from all chromosomes in the tophat .fai file to find the difference
	my @isect_ref = ();
	my @diff_ref = ();
	my %count = ();

	foreach my $element (@tophat_all_seqs, @isect) {
      $count{$element}++;
	};
	foreach my $element (keys %count) {
       push @{ $count{$element} > 1 ? \@isect_ref : \@diff_ref }, $element;
	};

	# Return the formatted exclusion list for the tophat fusion parameter --fusion-ignore-chromosomes
	return join(",",@diff_ref);

}

sub filter_fusions {
	my $options = shift;

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my $sample = $options->{'sample'};
	my $post_outdir = File::Spec->catdir($options->{'tmp'}, 'tophatpostrun/tophatfusion_'.$sample);
	my $fusions_file = File::Spec->catfile($post_outdir, 'result.txt');
	die "Please run tophatfusion_post step prior to filter\n" unless(-d $post_outdir);
	die "Please run tophatfusion_post step prior to filter\n" unless(-e $fusions_file);

	if(-s $fusions_file == 0) {
		system("echo '##EOF##' > $fusions_file") && die "An error occurred: $!";
	}

	my $normals_file = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},'cgpRna',$options->{'normalfusionslist'});

	my $command = "$^X ";
	$command .= _which('filter_fusions.pl');
	$command .= sprintf $FUSIONS_FILTER, 	$fusions_file,
						$sample,
						$normals_file,
						$post_outdir;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	return 1;
}

sub process_rg_tags {
	my $options = shift;

	my $rg_line = $options->{'rgline'};
	$rg_line =~ s/'//g;
	my @rg_fields = split(/\\t/, $rg_line);

	for(my $i=1 ; $i < scalar @rg_fields ; $i++){
		my @tag = split(/:/, $rg_fields[$i]);
		$options->{"$tag[0]"} = $tag[1];
	}
	return 1;
}

sub process_tophat_params {
	my ($options, $fusion_mode) = @_;

	my $ini_file = $options->{'config'};
	my $cfg = new Config::IniFiles( -file => $ini_file ) or die "Could not open config file: $ini_file";

	die "The tophat default parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($TOPHAT_DEFAULTS_SECTION));

	# Get the TopHat sample specific parameters from the options hash
	my $threads = $TOPHAT_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $TOPHAT_MAX_CORES);
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'num-threads', $threads);
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'output-dir', $options->{'tmp'}.'/tophat_'.$options->{'sample'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'transcriptome-index', $options->{'transcriptomepath'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'bowtie1', 'Y') if($options->{'bowtieversion'} != 2);
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'library-type', $options->{'librarytype'}) if(defined $options->{'librarytype'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-id', $options->{'ID'}) if(defined $options->{'ID'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-sample', $options->{'SM'}) if(defined $options->{'SM'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-library', $options->{'LB'}) if(defined $options->{'LB'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-description', '"'.$options->{'DS'}.'"') if(defined $options->{'DS'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-platform-unit', $options->{'PU'}) if(defined $options->{'PU'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-center', $options->{'CN'}) if(defined $options->{'CN'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-date', $options->{'DT'}) if(defined $options->{'DT'});
	$cfg->setval($TOPHAT_DEFAULTS_SECTION, 'rg-platform', $options->{'PL'}) if(defined $options->{'PL'});

	my @tophat_command;
	my @tophat_defaults = $cfg->Parameters($TOPHAT_DEFAULTS_SECTION);

	for my $key(@tophat_defaults){
		if($cfg->val($TOPHAT_DEFAULTS_SECTION, $key) eq 'Y'){
			push @tophat_command, "--".$key;
		}
		elsif($cfg->val($TOPHAT_DEFAULTS_SECTION, $key) ne '' && $cfg->val($TOPHAT_DEFAULTS_SECTION, $key) ne 'N'){
			push @tophat_command, "--".$key." ".$cfg->val($TOPHAT_DEFAULTS_SECTION, $key);
		}
	}

	if(defined $fusion_mode) {
		die "The tophat-fusion default parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($TOPHAT_FUSION_SECTION));
		$cfg->setval($TOPHAT_FUSION_SECTION, 'fusion-ignore-chromosomes', $options->{'exclude'}) if(defined $options->{'exclude'});

		my @tophat_fusion_defaults = $cfg->Parameters($TOPHAT_FUSION_SECTION);
		for my $key(@tophat_fusion_defaults){
			if($cfg->val($TOPHAT_FUSION_SECTION, $key) eq 'Y'){
				push @tophat_command, "--".$key;
			}
			elsif($cfg->val($TOPHAT_FUSION_SECTION, $key) ne '' && $cfg->val($TOPHAT_FUSION_SECTION, $key) ne 'N'){
				push @tophat_command, "--".$key." ".$cfg->val($TOPHAT_FUSION_SECTION, $key);
			}
		}

	}

	return join(" ", @tophat_command);
}

sub process_tophatpost_params {
	my $options = shift;

	my $ini_file = $options->{'config'};
	my $cfg = new Config::IniFiles( -file => $ini_file ) or die "Could not open config file: $ini_file";

	die "The tophat post parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($TOPHAT_POST_SECTION));

	my $threads = $TOPHAT_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $TOPHAT_MAX_CORES);
	$cfg->setval($TOPHAT_POST_SECTION, 'num-threads', $threads);
	$cfg->setval($TOPHAT_POST_SECTION, 'output-dir', './tophatfusion_'.$options->{'sample'});

	my @tophatpost_command;
	my @tophatpost_defaults = $cfg->Parameters($TOPHAT_POST_SECTION);
	for my $key(@tophatpost_defaults){
		if($cfg->val($TOPHAT_POST_SECTION, $key) eq 'Y'){
			push @tophatpost_command, "--".$key;
		}
		elsif($cfg->val($TOPHAT_POST_SECTION, $key) ne '' && $cfg->val($TOPHAT_POST_SECTION, $key) ne 'N'){
			push @tophatpost_command, "--".$key." ".$cfg->val($TOPHAT_POST_SECTION, $key);
		}
	}

	return join(" ", @tophatpost_command);
}

sub prog_version {
	my $options = shift;
	my $bowtie_path = $options->{'bowtiepath'};
	my $tophat_path = $options->{'tophatpath'};

	if(! defined $options->{'bowtiepath'} || $options->{'bowtiepath'} eq ''){
		my $release = '';
		$release = 2 if($options->{'bowtieversion'} == 2);
		$bowtie_path = _which('bowtie'.$release);
		$options->{'bowtiepath'} = $bowtie_path;
	}
	if(! defined $options->{'tophatpath'} || $options->{'tophatpath'} eq ''){
		$tophat_path = _which('tophat');
		$options->{'tophatpath'} = $tophat_path;
	}
	my $bowtie_version;
	my $tophat_version;
	{
		no autodie qw(system);
		my ($stdout, $stderr, $exit) = capture{ system("$bowtie_path --version"); };
		($bowtie_version) = $stdout =~ /version ([[:digit:]\.]+)/m;
		($stdout, $stderr, $exit) = capture{ system("$tophat_path --version"); };
		($tophat_version) = $stdout =~ /v([[:digit:]\.]+)/m;
	}
	return ($bowtie_version, $tophat_version);
}

sub split_setup {
	my $options = shift;

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my $sample = $options->{'sample'};
	my $tophat_fusion_dir = File::Spec->catdir($options->{'tmp'}, 'tophat_'.$sample);
	die "Please run the tophatfusion step prior to split\n" unless(-d $tophat_fusion_dir);
	my $fusions_file = File::Spec->catfile($tophat_fusion_dir,'fusions.out');
	die "Some files are missing from the tophatfusion process, please check the tophat fusion output directory and run the process again if necessary.\n" unless( -e $fusions_file && -e File::Spec->catfile($tophat_fusion_dir,'accepted_hits.bam') && -e File::Spec->catfile($tophat_fusion_dir,'unmapped.bam'));
	PCAP::Cli::file_for_reading('fusions.out', $fusions_file);

	# Create the TopHat Post run directory
	my $post_rundir = File::Spec->catdir($options->{'tmp'}, 'tophatpostrun');
	make_path($post_rundir) unless(-d $post_rundir);

	my $refgene = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},'tophat',$options->{'refgene'});
	my $ensgene = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},'tophat',$options->{'genebuild'},$options->{'ensgene'});
	my $blast = File::Spec->catdir($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},'tophat',$options->{'blastdb'});
	symlink($refgene, $post_rundir.'/refGene.txt') unless(-l File::Spec->catfile($post_rundir,'refGene.txt'));
	symlink($ensgene, $post_rundir.'/ensGene.txt') unless(-l File::Spec->catfile($post_rundir,'ensGene.txt'));
	symlink($blast, $post_rundir.'/blast') unless(-l $post_rundir.'/blast');

	my @commands;
	# Reformat the fusions.out file created by tophatfusion to UCSC chromosome names if chr is not the chromosome name prefix already
	my $fusion_line = capture_stdout { system('head', '-n', 1, $fusions_file ); };
	if($fusion_line !~ /^chr/){
		push @commands, sprintf q{sed -i.bak -E 's/^([0-9|X|Y]+)-([0-9|X|Y]+.*)/chr\\1-chr\\2/' %s}, $fusions_file;
	}

	# Split the fusions file into sub-files
	push @commands, sprintf q{split --numeric-suffixes --verbose --lines=%s %s %s/fusions}, $FUSIONS_SPLIT, $fusions_file, $tophat_fusion_dir;

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);

	# Create the sample sub-directories in the Tophat post run directory based on the number of split fusion files
	my $split_count = 0;
	my $split_dir;
	my $file;
	my $bamlink;
	my $abs_tophat_rundir = File::Spec->rel2abs( $tophat_fusion_dir );
	opendir(my $dh, $abs_tophat_rundir);
	while(readdir $dh) {
		if($_ =~ m/^fusions[[:digit:]]+$/){
			$split_count++;
			$file = $_;
			$split_dir = File::Spec->catdir($post_rundir, 'tophat_'.$sample.'.'.$split_count);
			make_path($split_dir) unless(-d $split_dir);
			move($abs_tophat_rundir."/".$_ , $split_dir.'/fusions.out') or die "$!";
			symlink($abs_tophat_rundir.'/accepted_hits.bam', $split_dir.'/accepted_hits.bam') unless(-l File::Spec->catfile($split_dir,'accepted_hits.bam'));
			symlink($abs_tophat_rundir.'/junctions.bed', $split_dir.'/junctions.bed') unless(-l File::Spec->catfile($split_dir,'junctions.bed'));
		}
	}

	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	return 1;
}

sub tophat_fusion {
	my $options = shift;

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	my $sample = $options->{'sample'};

	# Check the input data
	my $input_meta = $options->{'meta_set'};

	# Get the RG header information to format parameters --rg-id and --rg-sample
	my $first_file = $input_meta->[0];
	my $rg_line;

	if($first_file->fastq) {
		$rg_line = q{'}.$first_file->rg_header(q{\t}).q{'};
	}
	else {
		($rg_line, undef) = PCAP::Bam::rg_line_for_output($first_file->in, $sample);
		$rg_line = q{'}.$rg_line.q{'};
	}
	$options->{'rgline'} = $rg_line;

	process_rg_tags($options);
	my $tophat_params = process_tophat_params($options, 1);

	my @input1;
	my @input2;
	for my $input(@{$input_meta}) {
		if($input->fastq) {
			# Paired fastq input
			if($input->paired_fq) {
				push @input1, $input->in.'_1.'.$input->fastq;
				push @input2, $input->in.'_2.'.$input->fastq;
			}
			# Interleaved fastq input
			else{
				push @input1, $input->in.'.'.$input->fastq;
			}
		}
		# BAM input. Paired fastq.gz files should be present in the tmp/input directory from the bamtofastq process.
		else{
			my $inputdir = File::Spec->catdir($tmp, 'input');
			my $abs_inputdir = File::Spec->rel2abs( $inputdir );
			my $rg = $input->rg;
			die "Please run the bamtofastq step prior to tophatfusion\n" unless(-e File::Spec->catfile($abs_inputdir, $sample.'.'.$rg.'_1.fastq.gz'));
			push @input1, File::Spec->catfile($abs_inputdir, $sample.'.'.$rg.'_1.fastq.gz');
			push @input2, File::Spec->catfile($abs_inputdir, $sample.'.'.$rg.'_2.fastq.gz');
		}
	}

	# Update the environment to use the correct version of bowtie
	my $bwtpath = dirname($options->{'bowtiepath'} );
	$ENV{PATH} = "$bwtpath:$ENV{PATH}" if($ENV{'PATH'} !~ /$bwtpath/);

	my $ref_index_stem = $options->{'referencepath'};
	my $tophat_path = $options->{'tophatpath'};
	if(! defined $tophat_path || $tophat_path eq ''){
	  $tophat_path = _which('tophat');
	}

	my $command = $tophat_path." ".$tophat_params." ".$ref_index_stem." ".join(",",@input1)." ".join(",",@input2);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);

	return 1;
}

sub tophatfusion_post {
	my $options = shift;

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);

	# Check that the tophat fusion post run directory exists and has been setup correctly
	my $post_rundir = File::Spec->catdir($options->{'tmp'}, 'tophatpostrun');
	die "Please run tophatfusion and split steps prior to tophatfusion_post\n" unless(-d $post_rundir);
	die "Some Tophat-fusion-post setup files are missing, please run tophatfusion and split steps prior to tophatfusion_post\n" unless( -l $post_rundir.'/ensGene.txt' && -l $post_rundir.'/refGene.txt');

	my $tophatpost_params = process_tophatpost_params($options);
	my $tophatpost = $options->{'tophatpath'};
	if(! defined $tophatpost || $tophatpost eq ''){
	  $tophatpost = _which('tophat-fusion-post');
	}
	my $tophatpostindex = $options->{'tophatpostpath'};

	my $runcommand = $tophatpost.'-fusion-post'." ".$tophatpost_params." ".$tophatpostindex;

	# Get the full path of the tophat post run directory as need to cd to that location immediately prior to running tophat post
	my $abs_path = File::Spec->rel2abs( $post_rundir );

	my $command = "cd $abs_path; $runcommand";

	# Ensure the correct version of bowtie is on the path along with blastn
	my $bwtpath = dirname($options->{'bowtiepath'} );
	my $blastnpath = $options->{'blastn'};

	if(! defined $options->{'blastn'} || $options->{'blastn'} eq ''){
	  $blastnpath = _which('blastn');
	  $options->{'blastn'} = $blastnpath;
	}

	$ENV{PATH} = "$bwtpath:$ENV{PATH}" if($ENV{'PATH'} !~ /$bwtpath/);
	$ENV{PATH} = "$blastnpath:$ENV{PATH}" if($ENV{'PATH'} !~ /$blastnpath/);
	_which('bowtie');

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);

	#If the output is empty, ensure that it passes further checks by adding ##EOF## to the file.
	my $output_post = $post_rundir.'/result.txt';
	if(!-e $output_post || -s $output_post == 0){
		system("echo '##EOF##' > $output_post") && die "An error occurred: $!";
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

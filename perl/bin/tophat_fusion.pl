#!/usr/bin/perl
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
use FindBin;
use lib "$FindBin::Bin/../lib";

use strict;
use warnings FATAL => 'all';
use autodie qw(:all);
use English qw( -no_match_vars );

use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use List::Util qw(first);
use Const::Fast qw(const);
use File::Copy;
use Config::IniFiles;
use version;
use Cwd;

use PCAP::Cli;
use Sanger::CGP::Tophat::Implement;

my $ini_file = '/nfs/users/nfs_a/am26/git/rnaseqEvaluation/perl/config/tophat.ini'; # default config.ini file path
const my @REQUIRED_PARAMS => qw(outdir sample);
const my @VALID_PROCESS => qw(bamtofastq tophatfusion split tophatpost filter);
const my %INDEX_FACTOR => (	'bamtofastq' => -1,
				'tophatfusion' => 1,
				'split' => 1,
				'tophatpost' => 1,
				'filter' => 1);

{
	my $options = setup();

	# bam_to_fastq will only be called if bam input is detected in the setup subroutine. The process is nulti-threaded so that multiple BAMs can be converted to fastq in parallel.
	if(exists $options->{'bam'} && (!exists $options->{'process'} || $options->{'process'} eq 'bamtofastq')){
		my $threads = PCAP::Threaded->new($options->{'threads'});
		&PCAP::Threaded::disable_out_err if(exists $options->{'index'});
		$threads->add_function('bamtofastq', \&Sanger::CGP::Tophat::Implement::bam_to_fastq);
		$threads->run($options->{'max_split'}, 'bamtofastq', $options);
	}
	
	Sanger::CGP::Tophat::Implement::tophat_fusion($options) if(!exists $options->{'process'} || $options->{'process'} eq 'tophatfusion');
	Sanger::CGP::Tophat::Implement::split_setup($options) if(!exists $options->{'process'} || $options->{'process'} eq 'split');
	Sanger::CGP::Tophat::Implement::tophatfusion_post($options)if(!exists $options->{'process'} || $options->{'process'} eq 'tophatpost');
	
	if(!exists $options->{'process'} || $options->{'process'} eq 'filter') {
		Sanger::CGP::Tophat::Implement::filter_fusions($options);
		cleanup($options);
	}
}

sub cleanup {
	my $options = shift;
	my $tmpdir = $options->{'tmp'};
	my $sample = $options->{'sample'};
	my $fusion_outdir = File::Spec->catdir($tmpdir, "tophat_$sample");
	my $post_outdir = File::Spec->catdir($options->{'tmp'}, 'tophatpostrun/tophatfusion_'.$sample);
	system("cp $post_outdir/result.html $options->{outdir}/$sample.tophatfusion.html");
	move(File::Spec->catfile($fusion_outdir, 'accepted_hits.bam'), $options->{outdir}) || die $!;
	move(File::Spec->catfile($fusion_outdir, 'unmapped.bam'), $options->{outdir}) || die $!;
	move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
	#remove_tree $tmpdir if(-e $tmpdir);
	return 0;
}

sub setup {
	my %opts;
	pod2usage(-msg => "\nERROR: Options must be defined.\n", -verbose => 1, -output => \*STDERR) if(scalar @ARGV == 0);
	$opts{'cmd'} = join " ", $0, @ARGV;
	
	GetOptions( 	'h|help' => \$opts{'h'},
			'm|man' => \$opts{'m'},
			'v|version' => \$opts{'version'},
			'o|outdir=s' => \$opts{'outdir'},
			's|sample=s' => \$opts{'sample'},
			'sp|species=s' => \$opts{'species'},
			'l|librarytype=s' => \$opts{'librarytype'},
			'rb|refbuild=s' => \$opts{'referencebuild'},
			'gb|genebuild=i' => \$opts{'genebuild'},
			'b|bowtie=i' => \$opts{'bowtieversion'},
			'r|refdataloc=s' => \$opts{'refdataloc'},
			'ri|refindex=s' => \$opts{'referenceindex'},
			'ti|transindex=s' => \$opts{'transcriptomeindex'},
			'ub|ucscbuild=s' => \$opts{'tophatpostbuild'},
			'ui|ucscindex=s' => \$opts{'tophatpostindex'},
			'n|normals=s' => \$opts{'normalfusionslist'},
			't|threads=i' => \$opts{'threads'},
			'p|process=s' => \$opts{'process'},
			'i|index=i' => \$opts{'index'},							
			'c|config=s' => \$opts{'config'},
	) or pod2usage(1);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
	pod2usage(-verbose => 2) if(defined $opts{'m'});

	# Read in the config.ini file
	$ini_file = $opts{'config'} if(defined $opts{'config'});
	die "No config file has been specified." if($ini_file eq '');
	my $cfg = new Config::IniFiles( -file => $ini_file ) or die "Could not open config file: $ini_file";
	$opts{'config'} = $ini_file;
	
	# Populate the options hash with values from the config file
	$opts{'refdataloc'} = $cfg->val('tophat-config','referenceloc') unless(defined $opts{'refdataloc'});
	$opts{'referencebuild'} = $cfg->val('tophat-config','referencebuild') unless(defined $opts{'referencebuild'});
	$opts{'genebuild'} = $cfg->val('tophat-config','genebuild') unless(defined $opts{'genebuild'});
	$opts{'bowtieversion'} = $cfg->val('tophat-config','bowtieversion') unless(defined $opts{'bowtieversion'});
	$opts{'referenceindex'} = $cfg->val('tophat-config','referenceindex') unless(defined $opts{'referenceindex'});
	$opts{'transcriptomeindex'} = $cfg->val('tophat-config','transcriptomeindex') unless(defined $opts{'transcriptomeindex'});
	$opts{'tophatpostbuild'} = $cfg->val('tophat-config','tophatpostbuild') unless(defined $opts{'tophatpostbuild'});
	$opts{'tophatpostindex'} = $cfg->val('tophat-config','tophatpostindex') unless(defined $opts{'tophatpostindex'});
	$opts{'ensgene'} = $cfg->val('tophat-config','ensgene') unless(defined $opts{'ensgene'});
	$opts{'refgene'} = $cfg->val('tophat-config','refgene') unless(defined $opts{'refgene'});
	$opts{'blastdb'} = $cfg->val('tophat-config','blastdb') unless(defined $opts{'blastdb'});
	$opts{'blastn'} = $cfg->val('tophat-config','blastn') unless(defined $opts{'blastn'});
	$opts{'normalfusionslist'} = $cfg->val('tophat-config','normalfusionslist') unless(defined $opts{'normalfusionslist'});
	$opts{'tophatpath'} = $cfg->val('tophat-config','tophatpath') unless(defined $opts{'tophatpath'});
	$opts{'species'} = $cfg->val('tophat-config','species') unless(defined $opts{'species'});
	$opts{'librarytype'} = $cfg->val('tophat-config','librarytype') unless(defined $opts{'librarytype'});

	# In case the bowtie version has been missed from the config.ini file and command line parameters, set it to the default value of 1
	$opts{'bowtieversion'} = 1 if(! defined $opts{'bowtieversion'} || $opts{'bowtieversion'} eq '');
	# If user entered, ensure the version of bowtie is 1 or 2
	PCAP::Cli::valid_index_by_factor('bowtieversion', $opts{'bowtieversion'}, 2, 1);
	
	# Lookup the bowtie program path to use based on the version of bowtie requested. The recommended default for TopHat Fusion is to use bowtie1
	if(defined $opts{'bowtieversion'} && $opts{'bowtieversion'} == 2){
		$opts{'bowtiepath'} = $cfg->val('tophat-config','bowtie2path');
	} 
	else {
		$opts{'bowtiepath'} = $cfg->val('tophat-config','bowtie1path');
	}
	
	my ($bowtie_version, $tophat_version) = Sanger::CGP::Tophat::Implement::prog_version(\%opts);

	# Print version information for this program as well as the bowtie and tophat algorithms being used
	if($opts{'version'}) {
		print 'CGP tophat_fusion.pl version: ',Sanger::CGP::Tophat::Implement->VERSION,"\n";
		print 'CGP bowtie version: ', $bowtie_version,"\n";
		print 'CGP tophat version: ', $tophat_version,"\n";
		exit 0;
	}
	die "bowtie can only be used with version 0.12.9+, the version found in path is: $bowtie_version\n" unless(version->parse($bowtie_version) >= version->parse('0.12.9'));
	die "tophat can only be used with version 2+, the version found in path is: $tophat_version\n" unless(version->parse($tophat_version) >= version->parse('2.0.0'));

	# Ensure the bowtieversion parameter is compatible with the actual version of the installed bowtie software
	if(version->parse($bowtie_version) >= version->parse('2.0.0') && $opts{'bowtieversion'} != 2){
		$opts{'bowtieversion'} = 2;
	}
		
	for(@REQUIRED_PARAMS) {
		pod2usage(-msg => "\nERROR: $_ is a required argument.\n", -verbose => 1, -output => \*STDERR) unless(defined $opts{$_});
	}
	
	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
	my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpTophatFusion');
	make_path($tmpdir) unless(-d $tmpdir);
	$opts{'tmp'} = $tmpdir;
	my $progress = File::Spec->catdir($tmpdir, 'progress');
	make_path($progress) unless(-d $progress);
	my $logs = File::Spec->catdir($tmpdir, 'logs');
	make_path($logs) unless(-d $logs);
	
	# Check the input is fastq (paired or interleaved) or BAM and that a mixture of these file types hasn't been entered
	$opts{'raw_files'} = \@ARGV;
	Sanger::CGP::Tophat::Implement::check_input(\%opts);
	
	if(exists $opts{'bam'}){
		my $input = File::Spec->catdir($tmpdir, 'input');
		make_path($input) unless(-d $input);
	}
	
	# Determine reference chromosomes to exclude from tophat mapping based on a comparison of the tophat reference and tophat post reference
	my $exclude_refs = Sanger::CGP::Tophat::Implement::check_ref_seqs(\%opts);
	$opts{'exclude'} = $exclude_refs;
	
	delete $opts{'process'} unless(defined $opts{'process'});
	delete $opts{'index'} unless(defined $opts{'index'});
	delete $opts{'config'} unless(defined $opts{'config'});

 	# Apply defaults
	$opts{'threads'} = 1 unless(defined $opts{'threads'});
	
	if(exists $opts{'process'}){
		PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
		my $max_index = $INDEX_FACTOR{$opts{'process'}};
		
		if($opts{'process'} eq 'bamtofastq'){
			die "Process bamtofastq is only valid for BAM input files\n" if(!exists $opts{'bam'});
			$max_index = $opts{'max_split'};
		}
		
		if(exists $opts{'index'}) {
			PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
			PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max_index, 1);
			$opts{'max_split'} = $opts{'index'};
		}
	}
	elsif(exists $opts{'index'}) {
	die "ERROR: -index cannot be defined without -process\n";
}

	return \%opts;
}

__END__

=head1 tophat_fusion.pl

Cancer Genome Project implementation of the TopHat Fusion RNA-Seq algorithm

=head1 SYNOPSIS

tophat_fusion.pl [options] [file(s)...]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name

  Optional
    -bowtie   		-b  	1 or 2. Whether to use bowtie1 or bowtie2 for the fusion search [1]
    -librarytype	-l  	Library type for the sample. Options for tophat are fr-unstranded, fr-firststrand or fr-secondstrand. [fr-unstranded]
    -threads   		-t  	Number of cores to use. [1]
    -config   		-c  	Path to config.ini file. Defaults for the reference and transcriptome related parameters are provided in the config.ini file.
    -refbuild 		-rb 	Reference assembly version. Can be UCSC or Ensembl format e.g. GRCh38 or hg38 [GRCh38] 
    -genebuild 		-gb 	Gene build version. This needs to be consistent with the reference build in terms of the version and chromosome name style [77]
    -refindex   	-ri 	Stem name of the bowtie index files for the reference which need to be compatible with the bowtie version [GRCh38.genome]
    -transindex		-ti 	Stem name of the bowtie index files for the transcriptome which need to be compatible with the bowtie version [GRCh38.77]
    -ucscbuild   	-ub 	Tophat fusion post requires a reference build in UCSC format which must be equivalent to the refbuild version specified e.g. if refbuild = GRCh37 ucscbuild should be hg19 [hg38]
    -ucscindex		-ui 	Stem name of the bowtie index files for the transcriptome in ucsc format which should be compatible with the bowtie version and ucsc build [hg38.genome]
    -species  		-sp 	Species [human]

  Targeted processing (further detail under OPTIONS):
    -process   		-p   	Only process this step then exit
    -index    		-i   	Only valid for process bamtofastq - 1..<num_input_bams>

  Other:
    -help      		-h   	Brief help message.
    -man       		-m   	Full documentation.
    -version   		-v   	Version

  File list can be full file names or wildcard, e.g.
    tophat_fusion.pl -t 16 -o myout -refbuild GRCh38 -genebuild 77 -s sample input/*.bam

  Run with '-m' for possible input file types and details on bowtie versions/index files.

=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  bamtofastq
  tophatfusion
  split
  tophatpost
  filter

=back

=head2 INPUT FILE TYPES

There are several types of file that the script is able to process.

=over 8

=item f[ast]q

A standard uncompressed fastq file. Requires a pair of inputs with standard suffix of '_1' and '_2'
immediately prior to '.f[ast]q' or an interleaved f[ast]q file where read 1 and 2 are adjacent
in the file.

=item f[ast]q.gz

As *.f[ast]q but compressed with gzip.

=item bam

A list of single lane BAM files, RG line is transfered to aligned files.

=back

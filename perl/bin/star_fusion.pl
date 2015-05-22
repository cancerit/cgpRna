#!/usr/bin/perl
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
use Sanger::CGP::Star::Implement;

use Data::Dumper;

my $ini_file = "$FindBin::Bin/../config/star.ini"; # default config.ini file path
const my @REQUIRED_PARAMS => qw(outdir sample);
const my @VALID_PROCESS => qw(prepare star starfusion filter);
const my %INDEX_FACTOR => (	'prepare' => -1,
				'star' => 1,
				'starfusion' => 1,
				'filter' => 1);

{
	my $options = setup();
	
print Dumper(\$options);
	
	if(!exists $options->{'process'} || $options->{'process'} eq 'prepare'){
		# Process the input files.
		my $threads = PCAP::Threaded->new($options->{'threads'});
		&PCAP::Threaded::disable_out_err if(exists $options->{'index'});
		$threads->add_function('prepare', \&Sanger::CGP::Star::Implement::prepare);
		$threads->run($options->{'max_split'}, 'prepare', $options);
	}

	Sanger::CGP::Star::Implement::star_chimeric($options) if(!exists $options->{'process'} || $options->{'process'} eq 'star');
	Sanger::CGP::Star::Implement::star_fusion($options) if(!exists $options->{'process'} || $options->{'process'} eq 'starfusion');
	
	if(!exists $options->{'process'} || $options->{'process'} eq 'filter'){
		Sanger::CGP::Star::Implement::filter_fusions($options);
		cleanup($options);
	}
}

sub cleanup {
	my $options = shift;
	my $tmpdir = $options->{'tmp'};
	my $star_outdir = File::Spec->catdir($options->{'tmp'}, 'star');
	Sanger::CGP::Star::Implement::sam_to_bam($options);
	move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs_star')) || die $!;
	move(File::Spec->catfile($star_outdir, 'Aligned.sortedByCoord.out.bam'), $options->{'outdir'}) || die $!;
	remove_tree $tmpdir if(-e $tmpdir);
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
			'rb|refbuild=s' => \$opts{'referencebuild'},
			'gb|genebuild=i' => \$opts{'genebuild'},
			'r|refdataloc=s' => \$opts{'refdataloc'},
			'g|gtffile=s' => \$opts{'gtffilename'},
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
	$opts{'refdataloc'} = $cfg->val('star-config','referenceloc') unless(defined $opts{'refdataloc'});
	$opts{'referencebuild'} = $cfg->val('star-config','referencebuild') unless(defined $opts{'referencebuild'});
	$opts{'genebuild'} = $cfg->val('star-config','genebuild') unless(defined $opts{'genebuild'});
	$opts{'normalfusionslist'} = $cfg->val('star-config','normalfusionslist') unless(defined $opts{'normalfusionslist'});
	$opts{'gtffilename'} = $cfg->val('star-config','gtffilename') unless(defined $opts{'gtffilename'});
	$opts{'species'} = $cfg->val('star-config','species') unless(defined $opts{'species'});
	$opts{'starpath'} = $cfg->val('star-config','starpath');
	$opts{'starfusionpath'} = $cfg->val('star-config','starfusionpath');

	# Print version information for this program
	if($opts{'version'}) {
		print 'CGP star_fusion.pl version: ',Sanger::CGP::Star::Implement->VERSION,"\n";
		my $star_version = Sanger::CGP::Star::Implement::prog_version(\%opts);
		print "Star version: $star_version\n";
		exit 0;
	}
		
	for(@REQUIRED_PARAMS) {
		pod2usage(-msg => "\nERROR: $_ is a required argument.\n", -verbose => 1, -output => \*STDERR) unless(defined $opts{$_});
	}
	
	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
	my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmpStar');
	make_path($tmpdir) unless(-d $tmpdir);
	$opts{'tmp'} = $tmpdir;
	my $progress = File::Spec->catdir($tmpdir, 'progress');
	make_path($progress) unless(-d $progress);
	my $logs = File::Spec->catdir($tmpdir, 'logs');
	make_path($logs) unless(-d $logs);
	my $inputdir = File::Spec->catdir($tmpdir, 'input');
	make_path($inputdir) unless(-d $inputdir);
	my $stardir = File::Spec->catdir($tmpdir, 'star');
	make_path($stardir) unless(-d $stardir);
	
	# Check the input is fastq (paired only) or BAM and that a mixture of these file types hasn't been entered
	$opts{'raw_files'} = \@ARGV;
	Sanger::CGP::Star::Implement::check_input(\%opts);

	delete $opts{'process'} unless(defined $opts{'process'});
	delete $opts{'index'} unless(defined $opts{'index'});
	delete $opts{'config'} unless(defined $opts{'config'});

 	# Apply defaults
	$opts{'threads'} = 1 unless(defined $opts{'threads'});
	
	if(exists $opts{'process'}){
		PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
		my $max_index = $INDEX_FACTOR{$opts{'process'}};
		
		$max_index = $opts{'max_split'} if($opts{'process'} eq 'prepare');
		
		if(exists $opts{'index'}) {
			PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
			if($opts{'process'} eq 'prepare'){
				PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, $max_index, 1);
				$opts{'max_split'} = $opts{'index'};
			}
			else{
				die "ERROR: -index is not a valid parameter for process $opts{'process'}, please re-run excluding the -index parameter.\n";
			}
		}
	}
	elsif(exists $opts{'index'}) {
	die "ERROR: -index cannot be defined without -process\n";
}

	return \%opts;
}

__END__

=head1 star_fusion.pl

Cancer Genome Project implementation of the STAR RNA-Seq algorithm
https://bitbucket.org/dranew/defuse

=head1 SYNOPSIS

star_fusion.pl [options] [file(s)...]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name

  Optional
    -gtffile 		-g  	GTF annotation file name which should be compatible with the refbuild and gene build versions. It should reside under /refdataloc/species/refbuild/genebuild/ [Homo_sapiens.GRCh38.77.gtf]
    -normals  	  	-n  	File containing list of gene fusions detected in normal samples using STAR. It should reside under /refdataloc/species/refbuild/normal-fusions/ [star-normal-fusions-b38]
    -threads   		-t  	Number of cores to use. [1]
    -config   		-c  	Path to config.ini file. It contains defaults for; the reference and gene build versions, star software and default star and star-fusion parameters [<cgpRna-install-location>/perl/config/star.ini]
    -refbuild 		-rb 	Reference assembly version. Can be UCSC or Ensembl format e.g. GRCh38 or hg38 [GRCh38] 
    -genebuild 		-gb 	Gene build version. This needs to be consistent with the reference build in terms of the version and chromosome name style. Please use the build number only minus any prefixes such as e/ensembl [77]
    -refdataloc  	-r  	Parent directory of the reference data
    -species  		-sp 	Species [human]

  Targeted processing (further detail under OPTIONS):
    -process   		-p   	Only process this step then exit
    -index    		-i   	Only valid for process prepare - 1..<num_input_files>

  Other:
    -help      		-h   	Brief help message.
    -man       		-m   	Full documentation.
    -version   		-v   	Version

  File list can be full file names or wildcard, e.g.
    star_fusion.pl -t 16 -o myout -refbuild GRCh38 -genebuild 77 -s sample input/*.bam

  Run with '-m' for possible input file types.

=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  prepare
  star
  starfusion
  filter

=back

=head2 INPUT FILE TYPES

There are several types of file that the script is able to process.

=over 8

=item f[ast]q

A standard uncompressed fastq file. Requires a pair of inputs with standard suffix of '_1' and '_2'
immediately prior to '.f[ast]q'.

=item f[ast]q.gz

As *.f[ast]q but compressed with gzip.

=item bam

A list of single lane BAM files, RG line is transfered to aligned files.

=back

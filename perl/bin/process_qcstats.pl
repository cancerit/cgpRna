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
use Const::Fast qw(const);
use Pod::Usage qw(pod2usage);
use version;
use Cwd;

use PCAP::Cli;
use Sanger::CGP::RnaQC::Implement;

const my @REQUIRED_PARAMS => qw(indir outdir sample);

{
	my $options = setup();
	
	# Reformat the rRNA output file so that it can be appended to the bam_stats output
	Sanger::CGP::RnaQC::Implement::rrna_stats($options);
	Sanger::CGP::RnaQC::Implement::read_distribution_stats($options);
	Sanger::CGP::RnaQC::Implement::gene_coverage($options);

	
}

sub setup {
	my %opts;
	pod2usage(-msg => "\nERROR: Options must be defined.\n", -verbose => 1, -output => \*STDERR) if(scalar @ARGV == 0);
	$opts{'cmd'} = join " ", $0, @ARGV;
	
	GetOptions( 	'h|help' => \$opts{'h'},
			'm|man' => \$opts{'m'},
			'v|version' => \$opts{'version'},
			'i|indir=s' => \$opts{'indir'},
			'o|outdir=s' => \$opts{'outdir'},
			's|sample=s' => \$opts{'sample'},
			't|threads=i' => \$opts{'threads'},
			'p|process=s' => \$opts{'process'},
	) or pod2usage(1);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
	pod2usage(-verbose => 2) if(defined $opts{'m'});
	
	# Print version information for this program
	if($opts{'version'}) {
		print 'CGP process_qcstats.pl version: ',Sanger::CGP::RnaQC::Implement->VERSION,"\n";
		exit 0;
	}
	
	for(@REQUIRED_PARAMS) {
		pod2usage(-msg => "\nERROR: $_ is a required argument.\n", -verbose => 1, -output => \*STDERR) unless(defined $opts{$_});
	}
	
	my $inputdir = $opts{'indir'};
	die "The input directory $inputdir does not exist, exiting...\n" unless(-d $inputdir);
	
	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
	# Check all of the QC stats files have been created in the input directory
	Sanger::CGP::RnaQC::Implement::check_input(\%opts);
	
	delete $opts{'process'} unless(defined $opts{'process'});
	
 	# Apply defaults
	$opts{'threads'} = 1 unless(defined $opts{'threads'});
	
	return \%opts;
}

__END__
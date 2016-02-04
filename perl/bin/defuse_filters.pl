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
use strict;
use warnings;

use autodie qw(:all);
use English qw( -no_match_vars );
use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use PCAP::Cli;

# Columns in the deFuse output file that will be filtered.
const my $DEFUSE_SPLIT_CHAR => '\t';

const my $SPLITR_MIN_PVAL_COL => 7;
const my $BREAKSEQS_ESTISLANDS_PERCIDENT_COL => 14;
const my $CDNA_BREAKSEQS_PERCIDENT_COL => 15;
const my $EST_BREAKSEQS_PERCIDENT_COL => 17;
const my $GENOME_BREAKSEQS_PERCIDENT_COL => 38;
const my $SPAN_COVERAGE_MIN_COL => 66;

const my $SPLITR_MIN_PVAL_VAL =>  0.1; # splitr_min_pvalue - > 0.1
const my $BREAKSEQS_ESTISLANDS_PERCIDENT_VAL => 0.3; # breakseqs_estislands_percident - < 0.3
const my $CDNA_BREAKSEQS_PERCIDENT_VAL => 0.1; # cdna_breakseqs_percident - < 0.1
const my $EST_BREAKSEQS_PERCIDENT_VAL => 0.3; # est_breakseqs_percident - < 0.3
const my $GENOME_BREAKSEQS_PERCIDENT_VAL => 0.1; # genome_breakseqs_percident - < 0.1
const my $SPAN_COVERAGE_MIN_VAL => 0.6; # span_coverage_min - > 0.6

{
	my $options = setup();

  my $input = File::Spec->rel2abs($options->{'input'});
  my $sample = $options->{'sample'};
  my $outdir = $options->{'outdir'};
  my $output = File::Spec->catfile($outdir, "$sample.defuse-fusion.normals.ext.filtered.txt");
  
  open (my $ifh, $input) or die "Could not open file '$input' $!";
	open(my $ofh, '>', $output) or die "Could not open file '$output' $!";
  
  while (<$ifh>) {
		chomp;
		my $line = $_;
		if($line =~ m/^breakpoint_ref/){
		  print $ofh $line."\tcgp_defuse_filter\n";
		}
		else{
		  my @fields = split $DEFUSE_SPLIT_CHAR, $line;
		  if($fields[$SPLITR_MIN_PVAL_COL-1] > $SPLITR_MIN_PVAL_VAL && $fields[$BREAKSEQS_ESTISLANDS_PERCIDENT_COL-1] < $BREAKSEQS_ESTISLANDS_PERCIDENT_VAL && 
		     $fields[$CDNA_BREAKSEQS_PERCIDENT_COL-1] < $CDNA_BREAKSEQS_PERCIDENT_VAL && $fields[$EST_BREAKSEQS_PERCIDENT_COL-1] < $EST_BREAKSEQS_PERCIDENT_VAL && 
		     $fields[$GENOME_BREAKSEQS_PERCIDENT_COL-1] < $GENOME_BREAKSEQS_PERCIDENT_VAL && $fields[$SPAN_COVERAGE_MIN_COL-1] > $SPAN_COVERAGE_MIN_VAL) {
		   
		    print $ofh $line."\t1\n";		     
		  }
		  else{
		    print $ofh $line."\t0\n";
		  }
		} 
  }
  close($ifh);
  close($ofh);

}

sub setup {
	my %opts;
	pod2usage(-msg => "\nERROR: Options must be defined.\n", -verbose => 1, -output => \*STDERR) if(scalar @ARGV == 0);
	$opts{'cmd'} = join " ", $0, @ARGV;
	
	GetOptions( 	'h|help' => \$opts{'h'},
			'm|man' => \$opts{'m'},
			'i|input=s' => \$opts{'input'},
			'o|outdir=s' => \$opts{'outdir'},
			's|sample=s' => \$opts{'sample'},			
	) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
	pod2usage(-verbose => 2) if(defined $opts{'m'});
	
	PCAP::Cli::file_for_reading('input', $opts{'input'});
		
	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
	return \%opts;
}

__END__

=head1 defuse_fusions.pl

Adds a flag (called cgp_defuse_filter) to the raw defuse data based on validation carried out by Graham Bignell on the CTTV RNA-Seq cell lines data set. 
The flag can be used to filter the data in downstream analysis with the aim of reducing the number of false positive fusions called. Details of the filter thresholds can be found in the constants section at the top of the script.

=head1 SYNOPSIS

defuse_fusions.pl [options]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name
    -input    		-i   	deFuse input file containing fusions called by the cgpRna pipeline.

In the output file, a row with a 1 in the column cgp_defuse_filter means that this fusion has passed the set of filter thresholds whereas 0 means this fusion can potentially be filtered out.
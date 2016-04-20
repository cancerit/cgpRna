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
use File::Copy;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use FindBin qw($Bin);
use File::Which qw(which);
use PCAP::Cli;
use PCAP::Threaded;

# Columns in the deFuse output file that will be filtered.
const my $GENERATE_ID => q{ awk '{print NR"\t"$0}' %s > %s};
const my $GRASS_INPUT => q{ awk 'OFS="\t" {if(NR!=1){print $18":"$20":"$19","$23":"$25":"$24, $4, $1}}' %s > %s};
const my $MT => q{ sed -i 's/M:/MT:/g' %s };
const my $RUN_GRASS => q{ -genome_cache %s -show_biotype -file %s };
const my $GRASS_FLAG => q{ awk '{print $3"\t"$NF}' %s > %s };
const my $SORT1 => q{ sort -k1,1 %s > %s };
const my $SORT2 => q{ sort -k1,1 %s > %s };
const my $JOIN => q{ join -t $'\t' %s %s > %s };
const my $HEADER => q{ head -n 1 %s | awk '{print $0"\tgrass_flag"}' > %s };
const my $COLLATE => q{ cat %s %s | cut -f2- > %s };

{
	my $options = setup();
  run_grass($options);
	add_flag($options);
	write_output($options);
	cleanup($options);
}

sub cleanup {
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my $sample = $options->{'sample'};
  move(File::Spec->catfile($tmp, $sample.'.infuse.detected.fusions.grass.txt'), $options->{'outdir'}) || die $!;
  remove_tree $tmp if(-e $tmp);
  return 0;
}

sub add_flag {
  my $options = shift;
  
  my $tmp = $options->{'tmp'};
  my $grass_file = File::Spec->catfile($tmp, 'fusion_data_grass_input_ann');
  if(-e $grass_file){
    my $command1 = sprintf $GRASS_FLAG, $grass_file, File::Spec->catfile($tmp, 'grass_flag');
    my $command2 = sprintf $SORT1, File::Spec->catfile($tmp, 'fusion_data_id'), File::Spec->catfile($tmp, 'fusion_data_id_sorted');
    my $command3 = sprintf $SORT2, File::Spec->catfile($tmp, 'grass_flag'),File::Spec->catfile($tmp, 'grass_flag_sorted');
    my $command4 = sprintf $JOIN, File::Spec->catfile($tmp, 'fusion_data_id_sorted'), File::Spec->catfile($tmp, 'grass_flag_sorted'), File::Spec->catfile($tmp, 'joined_fusion_data');
    my @commands = ($command1,$command2,$command3,$command4);
    PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
  }
  
  return 0;
}

sub run_grass {
  my $options = shift;
  my $input = File::Spec->rel2abs($options->{'input'});
  my $sample = $options->{'sample'};
  my $tmp = $options->{'tmp'};
  
  return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
  
  my $command1 = sprintf $GENERATE_ID, $input, File::Spec->catfile($tmp, 'fusion_data_id');
  my $command2 = sprintf $GRASS_INPUT, File::Spec->catfile($tmp, 'fusion_data_id'), File::Spec->catfile($tmp, 'fusion_data_grass_input');
  my $command3 = sprintf $MT, File::Spec->catfile($tmp, 'fusion_data_grass_input');
  my $command4 = _which('grass.pl');
  $command4 .= sprintf $RUN_GRASS, $options->{'cache'}, File::Spec->catfile($tmp, 'fusion_data_grass_input');
  
  my @commands = ($command1,$command2,$command3, $command4);
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
  PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
  
	return 0;
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
			'c|cache=s' => \$opts{'cache'},
	) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
	pod2usage(-verbose => 2) if(defined $opts{'m'});
	
	PCAP::Cli::file_for_reading('input', $opts{'input'});
		
	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
	# Create working directory for storing intermediate files
  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmp_grass');
  make_path($tmpdir) unless(-d $tmpdir);
  my $logdir = File::Spec->catdir($tmpdir, 'logs');
  make_path($logdir) unless(-d $logdir);
  my $progressdir = File::Spec->catdir($tmpdir, 'progress');
  make_path($progressdir) unless(-d $progressdir);
  $opts{'tmp'} = $tmpdir;
	
	return \%opts;
}

sub write_output {  
  my $options = shift;
  my $tmp = $options->{'tmp'};
  my $sample = $options->{'sample'};
  
  my $command1 = sprintf $HEADER, File::Spec->catfile($tmp, 'fusion_data_id'), File::Spec->catfile($tmp, 'out_header');
  my $command2 = sprintf $COLLATE, File::Spec->catfile($tmp, 'out_header'), File::Spec->catfile($tmp, 'joined_fusion_data'), File::Spec->catfile($tmp, $sample.'.infuse.detected.fusions.grass.txt');
  my @commands = ($command1,$command2);
  PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);  
  
  return 0;
}

sub _which {
  my $prog = shift;
  my $l_bin = $Bin;
  my $path = File::Spec->catfile($l_bin, $prog);
  $path = which($prog) unless(-e $path);
  die "Failed to find $prog in path or local bin folder ($l_bin)\n\tPATH: $ENV{PATH}\n" unless(defined $path && -e $path);
  return $path;
}

__END__

=head1 addGrassFlag.pl

Adds Grass fusion prediction flag to fusion data generated by CGP InFuse pipeline. Grass flag values are described here: https://github.com/cancerit/grass/wiki/Fusion-Flag-Value-Descriptions

=head1 SYNOPSIS

addGrassFlag.pl [options]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name
    -cache   		-c   	VAGrENT cache file
    -input    		-i   	Fusion data text file generated by the CGP InFuse pipeline (compare_overlapping_fusions.pl script).
    

The Grass flag is appended to the end of each row of the input file.
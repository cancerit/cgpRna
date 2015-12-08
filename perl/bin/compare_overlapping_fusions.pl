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
BEGIN {
  use Cwd qw(abs_path);
  use File::Basename;
  unshift (@INC,dirname(abs_path($0)).'/../lib');
};

use strict;
use warnings;

use autodie qw(:all);
use English qw( -no_match_vars );
use File::Path qw(remove_tree make_path);
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use File::Which qw(which);
use PCAP::Cli;
use File::Copy;
use PCAP::Threaded;
use Sanger::CGP::CompareFusions::Implement;
use Sanger::CGP::CompareFusions::FusionAnnotation;

use Data::Dumper;

const my @REQUIRED_PARAMS => qw(outdir sample gtf);
const my @VALID_PROCESS => qw(createjunctionbed runbedpairtopair processoverlaps singletons queryvagrent annotatebed selectannotation collateannotation deduplicate output);
const my %INDEX_FACTOR => (	'createjunctionbed' => -1,
				'runbedpairtopair' => 1,
				'processoverlaps' => 1,
				'singletons' => 1,
				'queryvagrent' => 1,
				'annotatebed' => 1,
				'selectannotation' => 1,
				'collateannotation' => 1,
				'deduplicate' => 1,
				'output' => 1);				
{
  my $options = setup();
  
  my $threads = PCAP::Threaded->new($options->{'threads'});
  &PCAP::Threaded::disable_out_err if(exists $options->{'index'});
	
  $threads->add_function('createjunctionbed', \&Sanger::CGP::CompareFusions::Implement::create_junction_bedpe);
	
  $threads->run($options->{'num'}, 'createjunctionbed', $options) if(!exists $options->{'process'} || $options->{'process'} eq 'createjunctionbed');
  
  Sanger::CGP::CompareFusions::Implement::run_bed_pairtopair($options) if(!exists $options->{'process'} || $options->{'process'} eq 'runbedpairtopair');
  Sanger::CGP::CompareFusions::Implement::process_overlap_files($options) if(!exists $options->{'process'} || $options->{'process'} eq 'processoverlaps');
  Sanger::CGP::CompareFusions::Implement::process_singletons($options) if(!exists $options->{'process'} || $options->{'process'} eq 'singletons');
  Sanger::CGP::CompareFusions::Implement::query_vagrent($options) if(!exists $options->{'process'} || $options->{'process'} eq 'queryvagrent');
  if(-s File::Spec->catfile($options->{'tmp'}, $options->{'sample'}.".1.bed") || -s File::Spec->catfile($options->{'tmp'}, $options->{'sample'}.".2.bed")){
    Sanger::CGP::CompareFusions::Implement::annotate_bed($options) if(!exists $options->{'process'} || $options->{'process'} eq 'annotatebed');
    Sanger::CGP::CompareFusions::Implement::select_annotation($options) if(!exists $options->{'process'} || $options->{'process'} eq 'selectannotation');
    Sanger::CGP::CompareFusions::Implement::collate_annotation($options) if(!exists $options->{'process'} || $options->{'process'} eq 'collateannotation');
    Sanger::CGP::CompareFusions::Implement::deduplicate_fusions($options) if(!exists $options->{'process'} || $options->{'process'} eq 'deduplicate');
  }
  
  if(!exists $options->{'process'} || $options->{'process'} eq 'output') {
    Sanger::CGP::CompareFusions::Implement::generate_output($options);
    cleanup($options);
  } 
}

sub cleanup {
  my $options = shift;
  my $tmpdir = $options->{'tmp'};
  my $sample = $options->{'sample'};
  move(File::Spec->catfile($tmpdir, "$sample.detected.fusions.txt"), $options->{'outdir'}) || die $!;
  move(File::Spec->catdir($tmpdir, 'logs'), File::Spec->catdir($options->{'outdir'}, 'logs')) || die $!;
  remove_tree $tmpdir if(-e $tmpdir);
  return 0;
}

sub setup {
  my %opts;
  pod2usage(-msg => "\nERROR: Options must be defined.\n", -verbose => 1, -output => \*STDERR) if(scalar @ARGV == 0);
  $opts{'cmd'} = join " ", $0, @ARGV;
	
  GetOptions( 	'h|help' => \$opts{'h'},
		'm|man' => \$opts{'m'},
		'g|gtf=s' => \$opts{'gtf'},
		'o|outdir=s' => \$opts{'outdir'},
		's|sample=s' => \$opts{'sample'},
		't|threads=i' => \$opts{'threads'},
		'p|process=s' => \$opts{'process'},
		'i|index=i' => \$opts{'index'},
		'c|cache=s' => \$opts{'cache'},
  ) or pod2usage(2);

  pod2usage(-verbose => 1) if(defined $opts{'h'});
  pod2usage(-verbose => 2) if(defined $opts{'m'});
	
  for(@REQUIRED_PARAMS) {
    pod2usage(-msg => "\nERROR: $_ is a required argument.\n", -verbose => 1, -output => \*STDERR) unless(defined $opts{$_});
  }
	
  PCAP::Cli::file_for_reading('gtf', $opts{'gtf'});
		
  # Check the output directory exists and is writeable, create if not
  PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
  # Create working directory for storing intermediate files
  my $tmpdir = File::Spec->catdir($opts{'outdir'}, 'tmp_compareFusions');
  make_path($tmpdir) unless(-d $tmpdir);
  my $logdir = File::Spec->catdir($tmpdir, 'logs');
  make_path($logdir) unless(-d $logdir);
  my $progressdir = File::Spec->catdir($tmpdir, 'progress');
  make_path($progressdir) unless(-d $progressdir);
  $opts{'tmp'} = $tmpdir;
	
  # Count the number of input files being compared and check the format is valid i.e. it's been generated by tophat, star or defuse
  my $file_count = scalar @ARGV;
  die "Please enter either two or three input files to compare." if($file_count <2 || $file_count > 3);
  $opts{'num'} = $file_count;
  $opts{'input_files'} = \@ARGV;
	
  my $format;
  my $format_num;
  my %fusion_files;
  my $input;
  for (my $iter=1; $iter <= $file_count; $iter++) {
    $input = $ARGV[$iter-1];
    $format = Sanger::CGP::CompareFusions::Implement::check_input($input);
    if($format eq 'star'){
      $format_num = 1;
    }
    elsif($format eq 'tophat'){
      $format_num = 2;
    }
    else{
      $format_num = 3;
    }
    $fusion_files{$format_num}{'format'} = $format;
    $fusion_files{$format_num}{'name'} = $input;
  }
	
  $opts{'fusion_files'} = \%fusion_files;
  delete $opts{'process'} unless(defined $opts{'process'});
  delete $opts{'index'} unless(defined $opts{'index'});
	
  $opts{'threads'} = 1 unless(defined $opts{'threads'});

  if(exists $opts{'process'}) {
    PCAP::Cli::valid_process('process', $opts{'process'}, \@VALID_PROCESS);
    my $max_index = $INDEX_FACTOR{$opts{'process'}};
    if($max_index == -1) {
      $max_index = $opts{'num'};
    }
    $opts{'max_index'} = $max_index;
    if(exists $opts{'index'}) {
      PCAP::Cli::opt_requires_opts('index', \%opts, ['process']);
      PCAP::Cli::valid_index_by_factor('index', $opts{'index'}, \@ARGV, $max_index);
    }
  }
  elsif(exists $opts{'index'}) {
    die "ERROR: -index cannot be defined without -process\n";
  }

  return \%opts;
}

__END__

=head1 compare_overlapping_fusions.pl

Produces a report of overlapping fusions that have been called by star-fusion and deFuse.

=head1 SYNOPSIS

compare_overlapping_fusions.pl [options] [file(s)...]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name
    -gtf    		-g   	GTF file to use with bedtools to annotate each fusion breakpoint position.
    -cache    		-c   	VAGrENT cache file that should be the same reference and gene build as the GTF file being used e.g. GRCh38 e77.
    
  Optional:
    -threads    	-t   	Number of threads (cpus) to use [1].
    
  Targeted processing (further detail under OPTIONS):
    -process   		-p   	Only process this step then exit
    -index    		-i   	Valid for processes; createbed, annotatebed and createbedpe - 1..<num_input_files>
    
    Input files should be in the format generated by cgpRna pipelines; defuse.pl or star_fusion.pl
    
=head1 OPTIONS

=over 2

=item B<-process>

Available processes for this tool are:

  createjunctionbed
  runbedpairtopair
  queryvagrent
  annotatebed
  selectannotation
  collateannotation
  output

=back

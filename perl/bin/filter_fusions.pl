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

use Data::Dumper;

# Position of the columns in the tophat-fusion-post output file used to format fusion breakpoint references.
const my $TOPHAT_SPLIT_CHAR => '\t';
const my $TOPHAT_CHR1 => 3;
const my $TOPHAT_POS1 => 4;
const my $TOPHAT_CHR2 => 6;
const my $TOPHAT_POS2 => 7;
const my @TOPHAT_HEADER => qw(breakpoint_ref sample_name gene_1 chr_1 pos_1 gene_2 chr_2 pos_2 num_spanning_reads num_spanning_mate_pairs num_spanning_mates_2 score);

# Position of the columns in the deFuse output file used to format fusion breakpoint references.
const my $DEFUSE_SPLIT_CHAR => '\t';
const my $DEFUSE_CHR1 => 25;
const my $DEFUSE_POS1 => 38;
const my $DEFUSE_CHR2 => 26;
const my $DEFUSE_POS2 => 39;
const my $DEFUSE_HEADER_PATTERN => '^cluster_id';

# Position of the columns in the star-fusion output file used to format fusion breakpoint references.
const my $STAR_SPLIT_CHAR => '\t';
const my $STAR_SPLIT_CHAR2 => ':';
const my $STAR_CHR1 => 5;
const my $STAR_CHR2 => 8;
const my $STAR_HEADER_PATTERN => '^#fusion_name';

{
	my $options = setup();

	reformat($options);
	normals_filter($options);
	write_output($options);
	remove_tree $options->{'tmp'} if(-e $options->{'tmp'});
}

sub normals_filter {
	my $options = shift;
	my $tmp = $options-> {'tmp'};
	my $sample = $options-> {'sample'};
	my $fusions_file = File::Spec->catfile($tmp,"$sample.fusions");
	my $fusions_sorted = File::Spec->catfile($tmp,"$sample.fusions.sorted");
	PCAP::Cli::file_for_reading('fusions.reformat', $fusions_file.'.reformat');
	
	# Sort the fusions file prior to joining with the normals file
	system("sort -k1,1 $fusions_file.reformat > $fusions_file.sorted");
	
	my $normal_file = $options-> {'normals'};
	system("join -v 2 $normal_file $fusions_file.sorted > $fusions_file.filtered");
	
	return 1;
}

sub process_program_params {
	my $options = shift;
	
	my $program = lc($options->{'program'});
	if($program eq 'tophat'){
		$options->{'split-char'} = $TOPHAT_SPLIT_CHAR;
		$options->{'chr1'} = $TOPHAT_CHR1;
		$options->{'pos1'} = $TOPHAT_POS1;
		$options->{'chr2'} = $TOPHAT_CHR2;
		$options->{'pos2'} = $TOPHAT_POS2;
		$options->{'header'} = join("\t",@TOPHAT_HEADER);
	}
	elsif($program eq 'defuse'){
		$options->{'split-char'} = $DEFUSE_SPLIT_CHAR;
		$options->{'chr1'} = $DEFUSE_CHR1;
		$options->{'pos1'} = $DEFUSE_POS1;
		$options->{'chr2'} = $DEFUSE_CHR2;
		$options->{'pos2'} = $DEFUSE_POS2;
		$options->{'header-pattern'} = $DEFUSE_HEADER_PATTERN;
	}
	elsif($program eq 'star'){
		$options->{'split-char'} = $STAR_SPLIT_CHAR;
		$options->{'split-char2'} = $STAR_SPLIT_CHAR2;
		$options->{'chr1'} = $STAR_CHR1;
		$options->{'chr2'} = $STAR_CHR2;
		$options->{'header-pattern'} = $STAR_HEADER_PATTERN;
	}
	else{
		die "An invalid program name has been entered, parameter -p should be tophat, defuse or star.\n"
	}
	
	return 1;
}

sub reformat {
	my $options = shift;
	my $input = File::Spec->rel2abs($options->{'input'});
	my $tmp = $options->{'tmp'};
	my $sample = $options->{'sample'};
	my $output = File::Spec->catfile($tmp, "$sample.fusions.reformat");
	my $program = lc($options->{'program'});
	my $header_pattern = "###";

	open (my $ifh, $input) or die "Could not open file '$input' $!";
	open(my $ofh, '>', $output) or die "Could not open file '$output' $!";
	
	while (<$ifh>) {
		chomp;
		my $line = $_;
		
		$header_pattern = $options->{'header-pattern'} if(exists $options->{'header-pattern'});
		if($line =~ m/$header_pattern/){
			$options->{'header'} = "breakpoint_ref\t".$line;
		}
		else{
			$line =~ s/\tchr/\t/g;
			my @fields = split '\t', $line;
			my $fusion;
			
			if($program eq 'star'){
				my @break1 = split ':', $fields[$options->{'chr1'} - 1];
				my @break2 = split ':', $fields[$options->{'chr2'} - 1];
				$fusion = $break1[0].":".$break1[1]."-".$break2[0].":".$break2[1];
			}
			else{
				$fusion = $fields[$options->{'chr1'}-1].":".$fields[$options->{'pos1'}-1]."-".$fields[$options->{'chr2'}-1].":".$fields[$options->{'pos2'}-1];
			}
			
			print $ofh "$fusion\t$line\n";
		}
	}	
	close ($ifh);
	close ($ofh);
	
	return 1;
}

sub setup {
	my %opts;
	pod2usage(-msg => "\nERROR: Options must be defined.\n", -verbose => 1, -output => \*STDERR) if(scalar @ARGV == 0);
	$opts{'cmd'} = join " ", $0, @ARGV;
	
	GetOptions( 	'h|help' => \$opts{'h'},
			'm|man' => \$opts{'m'},
			'i|input=s' => \$opts{'input'},
			'o|outdir=s' => \$opts{'outdir'},
			'n|normals=s' => \$opts{'normals'},
			's|sample=s' => \$opts{'sample'},
			'p|program=s' => \$opts{'program'},						
	) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
	pod2usage(-verbose => 2) if(defined $opts{'m'});
	
	PCAP::Cli::file_for_reading('input', $opts{'input'});
	PCAP::Cli::file_for_reading('normals', $opts{'normals'});
	
	# Setup details specific to each fusion detection program output file
	process_program_params(\%opts);
		
	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});
	
	# Create working directory for storing intermediate files
	my $tmpdir = File::Spec->catdir($opts{'outdir'}, "tmp_$opts{'program'}Filter");
	make_path($tmpdir) unless(-d $tmpdir);
	$opts{'tmp'} = $tmpdir;

	return \%opts;
}

sub write_output {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	my $sample = $options-> {'sample'};
	my $program = $options-> {'program'};
	my $outdir = $options->{'outdir'};
	my $fusions_file = File::Spec->catfile($tmp,"$sample.fusions.filtered");
	my $output_file = File::Spec->catfile($outdir,"$sample.$program-fusion.normals.filtered.txt");
	PCAP::Cli::file_for_reading('filtered.fusions', $fusions_file);
	
	open (my $ifh, $fusions_file) or die "Could not open file $fusions_file $!";
	open(my $ofh, '>', $output_file) or die "Could not open file $output_file $!";
	print $ofh $options->{'header'}."\n";
	while (<$ifh>) {
		chomp;
		my $line = $_;
		$line =~ s/\s/\t/g;
		print $ofh $line."\n";
	}
	close ($ifh);
	close ($ofh);	
	
	return 1;
}

__END__

=head1 filter_fusions.pl

Reformats and filters the output file from star-fusion, tophatfusion-post or deFuse against a file containing gene fusions called in normal samples.

=head1 SYNOPSIS

filter_fusions.pl [options]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name
    -program   		-p   	Gene fusion detection program used, valid values are; tophat, star or defuse
    -input    		-i   	Input file of fusions generated by the program (e.g. for star; star-fusion.fusion_candidates.txt).
    -normals   		-n   	File containing list of gene fusions detected in normal samples using the relevant program (e.g. star-normal-fusions-b38)
                    		 Expected format is one column pre-sorted using the Unix sort command <chr1:pos1-chr2:pos2> e.g.
                    		  10:100000026-12:93371978
                    		  10:100000026-X:84180396
                    		  ...

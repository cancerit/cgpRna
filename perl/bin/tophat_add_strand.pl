#!/usr/bin/perl
##########LICENCE ##########
#Copyright (c) 2015-2017 Genome Research Ltd.
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
use Getopt::Long;
use File::Spec;
use Pod::Usage qw(pod2usage);
use Const::Fast qw(const);
use PCAP::Cli;

const my @TOPHAT_HEADER => qw(breakpoint_ref sample_name gene_1 chr_1 pos_1 gene_2 chr_2 pos_2 num_spanning_reads num_spanning_mate_pairs num_spanning_mates_2 score strand_1 strand_2);

{
	my $options = setup();

	reformat($options);
	add_strand($options);
	write_output($options);
}

sub add_strand {
	my $options = shift;
	my $outdir = $options->{'outdir'};
	my $sample = $options-> {'sample'};
	my $fusions_file = $options-> {'input'};
	my $strand_file = File::Spec->catfile($outdir,"$sample.tophat.strand");
	my $output = File::Spec->catfile($outdir,"$sample.tophat.fusions");

	PCAP::Cli::file_for_reading('filtered.fusions', $fusions_file);

	# Sort both files prior to joining
	system("sort -k1,1 $fusions_file > $output.sorted") && die "An error occurred: $!";
	system("sort -k1,1 $strand_file > $strand_file.sorted") && die "An error occurred: $!";
	system("join $output.sorted $strand_file.sorted > $output.strand") && die "An error occurred: $!";

	if(-s "$output.strand" == 0) {
		system("echo '##EOF##' > $output.strand") && die "An error occurred: $!";
	}

	return 1;
}

sub reformat {
	my $options = shift;
	my $potential_fusions = File::Spec->rel2abs($options->{'potfusions'});
	my $outdir = $options->{'outdir'};
	my $sample = $options->{'sample'};
	my $output = File::Spec->catfile($outdir, "$sample.tophat.strand");

	open (my $ifh, $potential_fusions) or die "Could not open file '$potential_fusions' $!";
	open(my $ofh, '>', $output) or die "Could not open file '$output' $!";

	while (<$ifh>) {
		chomp;
		my $line = $_;
		if($line =~ m/^$sample/){

			my @fields = split ' ', $line;
			$fields[1] =~ s/-|chr/ /g;
			my @strands = split '', $fields[4];
			my $strand_upd = join(" ",@strands);
			$strand_upd =~ tr/rf/-+/;
			$fields[4] = $strand_upd;
			$line = join(" ",@fields);
			my @fields_upd = split ' ', $line;
			my $pos1 = $fields_upd[3] + 1;
			my $pos2 = $fields_upd[4] + 1;
			print $ofh "$fields_upd[1]:$pos1-$fields_upd[2]:$pos2 $fields_upd[5] $fields_upd[6]\n";

		}
	}
	close ($ifh);
	close ($ofh);

	if(-s $output == 0) {
		system("echo '##EOF##' > $output") && die "An error occurred: $!";
	}

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
			's|sample=s' => \$opts{'sample'},
			'p|potfusions=s' => \$opts{'potfusions'},
	) or pod2usage(2);

	pod2usage(-verbose => 1) if(defined $opts{'h'});
	pod2usage(-verbose => 2) if(defined $opts{'m'});

	PCAP::Cli::file_for_reading('input', $opts{'input'});
	PCAP::Cli::file_for_reading('potential_fusions', $opts{'potfusions'});

	# Check the output directory exists and is writeable, create if not
	PCAP::Cli::out_dir_check('outdir', $opts{'outdir'});

	return \%opts;
}

sub write_output {
	my $options = shift;

	my $sample = $options-> {'sample'};
	my $outdir = $options->{'outdir'};
	my $fusions_file = File::Spec->catfile($outdir,"$sample.tophat.fusions".'.strand');
	my $output_file = File::Spec->catfile($outdir,"$sample.tophat-fusion.normals.filtered.strand.txt");
	PCAP::Cli::file_for_reading('filtered.fusions', $fusions_file);

	open (my $ifh, $fusions_file) or die "Could not open file $fusions_file $!";
	open(my $ofh, '>', $output_file) or die "Could not open file $output_file $!";
	print $ofh join("\t", @TOPHAT_HEADER)."\n";
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

=head1 tophat_add_strand.pl

Reformats and filters the tophat output from filter_fusions.pl, to add strand information. It requires the file potential_fusion.txt, produced by tophatfusion-post.

=head1 SYNOPSIS

tophat_add_strand.pl [options]

  Required parameters:
    -outdir    		-o   	Folder to output result to.
    -sample   		-s   	Sample name
    -potfusions   	-p   	location of the potential_fusion.txt file produced by tophatfusion-post
    -input   		-i   	File containing list of gene fusions that have been detected by tophatfusion-post and filtered by filter_fusions.pl
                    		 Expected format is tab separated with the first eight columns: breakpoint_ref	sample_name	gene_1	chr_1	pos_1	gene_2	chr_2	pos_2 e.g.
                    		  11:18247996-9:114324016	TEST.3	SAA2	11	18247996	ORM1	9	114324016
                    		  5:134927647-19:10396126	TEST.2	PCBD2	5	134927647	CDC37	19	10396126
                    		  ...

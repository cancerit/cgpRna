#!/usr/bin/perl

##########LICENCE ########
# Copyright (c) 2015-2016 Genome Research Ltd.
#
# Author: Cancer Genome Project <cgpit@sanger.ac.uk>
#
# This file is part of cgpRna.
#
# cgpRna is free software: you can redistribute it and/or modify it under
# the terms of the GNU Affero General Public License as published by the
# Free Software Foundation; either version 3 of the License, or (at your
# option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Affero
# General Public License for more details.
#
# You should have received a copy of the GNU Affero General Public
# License along with this program. If not, see
# <http://www.gnu.org/licenses/>.
#
# 1. The usage of a range of years within a copyright statement contained
# within this distribution should be interpreted as being equivalent to a
# list of years including the first and last year specified and all
# consecutive years between them. For example, a copyright statement that
# reads ‘Copyright (c) 2005, 2007- 2009, 2011-2012’ should be interpreted
# as being identical to a statement that reads ‘Copyright (c) 2005, 2007,
# 2008, 2009, 2011, 2012’ and a copyright statement that reads ‘Copyright
# (c) 2005-2012’ should be interpreted as being identical to a statement
# that reads ‘Copyright (c) 2005, 2006, 2007, 2008, 2009, 2010, 2011,
# 2012’."
##########LICENCE ##########

use strict;
use Bio::DB::HTS::Tabix;
use Data::Dumper;

## Required input
#SNP6 annotations using script in vagrent git annotateSV.pl
# Fusion gene data from RNA-seq pipeline


my ($snp_sv,$fusion_genes)=@ARGV;

print "\nUsage: perl compare_CN.pl SNP6_sv_annotated_data.out FusionGeneData.out\n" if( (!defined $snp_sv) || (!defined $fusion_genes) );

system("bgzip -f $snp_sv") if(! -s $snp_sv.'.gz');

my @header =qw/#chr	lStart	lEnd	chr	rStart	rEnd	name	score	lStrand	rStrand	microHomoLen	cancer_type	patientId	sampleId	lExon\/rExon	lTranscript\/rTranscript	lGeneStrand\/rGeneStrand	lAnnotation\/rAnnotation	lGene\/rGene	SpannedGenes	lPromoterL\/lPromoterR	rPromoterL\/rPromoterR	lEnhancerL\/lEnhancerR	rEnhancerL\/rEnhancerR/;

system("tabix -s 1 -b 2 -e 6 -p bed $snp_sv.gz") if(! -s $snp_sv.'.gz.tbi');

open my $fusion_fh, '<', $fusion_genes;

# create tabix file
my $tabix = Bio::DB::HTS::Tabix(filename => "$snp_sv.gz");

open my $matched_fh, '>', 'matched'.$snp_sv.'_and_'.$fusion_genes;

my $extend_cutoff=1000;

my $count=0;
while (<$fusion_fh>) {
	chomp;
	$count++;
	if($count == 1) {	print $matched_fh "Match[2-full,1-partial,0-none]\tExons\tbreakpoints\tcalledby\t".join("\t",@header)."\n"; next;}
	my $flag=0;
	# Note gene is transcript id
	my($exons,$breakpoint,$called,$chr1,$transcript1,$exon1_start,$exon1_end,$chr2,$transcript2,$exon2_start,$exon2_end)=(split "\t",$_)[0,1,2,3,6,8,9,10,13,15,16];
	my $fusion_line="$exons\t$breakpoint\t$called";
 # print "======>$exons,$breakpoint,$called,$chr1,$transcript1,$exon1_start,$exon1_end,$chr2,$transcript2,$exon2_start,$exon2_end \n";
  # For interchromosomal translocations extend the exon boundaries by $extend_cutoff paramater.
  if($chr1 ne $chr2) {
		($exon1_start,$exon1_end)=get_start_end($exon1_start,$exon1_end);
		my($results)=search_break_point_annotations($chr1,($exon1_start - $extend_cutoff),($exon1_end+$extend_cutoff),$tabix,$transcript1,"NAAA",$fusion_line);
		print_results($matched_fh,$results);
		($exon2_start,$exon2_end)=get_start_end($exon2_start,$exon2_end);
		$results=search_break_point_annotations($chr2,($exon2_start - $extend_cutoff),($exon2_end+$extend_cutoff),$tabix,$transcript2,"NAAA",$fusion_line);
	 	print_results($matched_fh,$results);
	}
	if($chr1 eq $chr2) {
		($exon1_start,$exon1_end)=get_start_end($exon1_start,$exon1_end);
		($exon2_start,$exon2_end)=get_start_end($exon2_start,$exon2_end);
		# final check to see start and end are in right direction
		($exon1_start,$exon2_end)=get_start_end($exon1_start,$exon2_end);
		my($results)=search_break_point_annotations($chr1,$exon1_start, $exon2_end,$tabix,$transcript1,$transcript2,$fusion_line);
		print_results($matched_fh,$results);
 }
 #print "*********$_\n $chr1,$exon1_start,$exon1_end ------- $chr2,$exon2_start,$exon2_end \n";
 #print $matched_fh "<<<<<<<< EXON Fusion Data >>>>>>>>> \n $line\n <<<<<<<< SNP ArrayCN Data >>>>>>>>>\nMatch[2-full,1-partial,0-none]\t$header_sv";


}
close($matched_fh);

sub print_results {
	my($fh,$results)=@_;
	if(defined $results ) {
		foreach my $line(@$results) {
			print $fh "$line\n";
		}
	}
}


sub get_start_end {
	my($start,$end)=@_;
	if($start > $end ) { my $tmp=$end; $end=$start; $start=$tmp;}
  return($start,$end);
}

sub search_break_point_annotations {
 my($chr,$start,$end,$tabix,$transcript1,$transcript2,$fusion_line)=@_;
	my $results=undef;
  my $iter = $tabix->query(sprintf '%s:%d-%d', $chr,$start,$end);
  while(my $record = $iter->next){
    my ($lgene,$rgene)=(split "/", (split "\t",$record)[12]);
    my($matched)=compare_transcripts($transcript1,$transcript2,$lgene,$rgene);
    push(@$results,"$matched\t$fusion_line\t$record\n") if($matched );
  }
return $results;
}

sub compare_transcripts {
	my($t1,$t2,$lg,$rg)=@_;
	my $matched=0;
		$matched++ if(($lg eq $t1) or ($lg eq $t2));
		$matched++ if(($rg eq $t1) or ($rg eq $t2));
	  return ($matched);
}


sub compare_gene {
	my($gene,$lg,$rg)=@_;
	my $matched=0;
# can use $matched++ if {first {$gene,$_}} @lg_array;
	my @lg_array=split(',',$lg);
	foreach my $g(@lg_array) {
		$matched++ if($g eq $gene);
 	}
	my @rg_array=split(',',$rg);
	foreach my $g(@rg_array) {
		$matched++ if($g eq $gene);
 	}
	return $matched;
}








package Sanger::CGP::CompareFusions::FusionAnnotation;
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

use strict;

use Sanger::CGP::CgpRna;
our $VERSION = Sanger::CGP::CgpRna->VERSION;

sub new {
  my ($class, %args) = @_;
  my $self = {};
  bless $self,$class;
  
  if ($args{-breakpoint}) { $self->breakpoint($args{-breakpoint}) }
  if ($args{-alt_breakpoint}) { $self->alt_breakpoint($args{-alt_breakpoint}) }
  if ($args{-chr1}) { $self->chr1($args{-chr1}) }
  if ($args{-strand1}) { $self->strand1($args{-strand1}) }
  if ($args{-pos1_start}) { $self->pos1_start($args{-pos1_start}) }
  if ($args{-pos1_end}) { $self->pos1_end($args{-pos1_end}) }
  if ($args{-feature1}) { $self->feature1($args{-feature1}) }
  if ($args{-feature1_start}) { $self->feature1_start($args{-feature1_start}) }
  if ($args{-feature1_end}) { $self->feature1_end($args{-feature1_end}) }
  if ($args{-exon1_num}) { $self->exon1_num($args{-exon1_num}) }
  if ($args{-exon1_id}) { $self->exon1_id($args{-exon1_id}) }
  if ($args{-gene1}) { $self->gene1($args{-gene1}) }
  if ($args{-gene1_id}) { $self->gene1_id($args{-gene1_id}) }
  if ($args{-gene1_start}) { $self->gene1_start($args{-gene1_start}) }
  if ($args{-gene1_end}) { $self->gene1_end($args{-gene1_end}) }
  if ($args{-transcript1_id}) { $self->transcript1_id($args{-transcript1_id}) }
  if ($args{-distance1}) { $self->distance1($args{-distance1}) }
  if ($args{-chr2}) { $self->chr2($args{-chr2}) }
  if ($args{-strand2}) { $self->strand2($args{-strand2}) }
  if ($args{-pos2_start}) { $self->pos2_start($args{-pos2_start}) }
  if ($args{-pos2_end}) { $self->pos2_end($args{-pos2_end}) }
  if ($args{-feature2}) { $self->feature2($args{-feature2}) }
  if ($args{-feature2_start}) { $self->feature2_start($args{-feature2_start}) }
  if ($args{-feature2_end}) { $self->feature2_end($args{-feature2_end}) }
  if ($args{-exon2_num}) { $self->exon2_num($args{-exon2_num}) }
  if ($args{-exon2_id}) { $self->exon2_id($args{-exon2_id}) }
  if ($args{-gene2}) { $self->gene2($args{-gene2}) }
  if ($args{-gene2_id}) { $self->gene2_id($args{-gene2_id}) }
  if ($args{-gene2_start}) { $self->gene2_start($args{-gene2_start}) }
  if ($args{-gene2_end}) { $self->gene2_end($args{-gene2_end}) }
  if ($args{-transcript2_id}) { $self->transcript2_id($args{-transcript2_id}) }
  if ($args{-distance2}) { $self->distance2($args{-distance2}) }
  
  return $self;
}

sub breakpoint {
  my $self = shift;
  $self->{breakpoint} = shift if @_;
  return($self->{breakpoint});
}

sub alt_breakpoint {
  my $self = shift;
  $self->{alt_breakpoint} = shift if @_;
  return($self->{alt_breakpoint});
}

sub chr1 {
  my $self = shift;
  $self->{chr1} = shift if @_;
  return($self->{chr1});
}

sub chr2 {
  my $self = shift;
  $self->{chr2} = shift if @_;
  return($self->{chr2});
}

sub pos1_start {
  my $self = shift;
  $self->{pos1_start} = shift if @_;
  return($self->{pos1_start});
}

sub pos2_start {
  my $self = shift;
  $self->{pos2_start} = shift if @_;
  return($self->{pos2_start});
}

sub pos1_end {
  my $self = shift;
  $self->{pos1_end} = shift if @_;
  return($self->{pos1_end});
}

sub pos2_end {
  my $self = shift;
  $self->{pos2_end} = shift if @_;
  return($self->{pos2_end});
}

sub strand1 {
  my $self = shift;
  $self->{strand1} = shift if @_;
  return($self->{strand1});
}

sub strand2 {
  my $self = shift;
  $self->{strand2} = shift if @_;
  return($self->{strand2});
}

sub feature1 {
  my $self = shift;
  $self->{feature1} = shift if @_;
  return($self->{feature1});
}

sub feature2 {
  my $self = shift;
  $self->{feature2} = shift if @_;
  return($self->{feature2});
}

sub feature1_start {
  my $self = shift;
  $self->{feature1_start} = shift if @_;
  return($self->{feature1_start});
}

sub feature2_start {
  my $self = shift;
  $self->{feature2_start} = shift if @_;
  return($self->{feature2_start});
}

sub feature1_end {
  my $self = shift;
  $self->{feature1_end} = shift if @_;
  return($self->{feature1_end});
}

sub feature2_end {
  my $self = shift;
  $self->{feature2_end} = shift if @_;
  return($self->{feature2_end});
}

sub gene1 {
  my $self = shift;
  $self->{gene1} = shift if @_;
  return($self->{gene1});
}

sub gene2 {
  my $self = shift;
  $self->{gene2} = shift if @_;
  return($self->{gene2});
}

sub gene1_id {
  my $self = shift;
  $self->{gene1_id} = shift if @_;
  return($self->{gene1_id});
}

sub gene2_id {
  my $self = shift;
  $self->{gene2_id} = shift if @_;
  return($self->{gene2_id});
}

sub gene1_start {
  my $self = shift;
  $self->{gene1_start} = shift if @_;
  return($self->{gene1_start});
}

sub gene2_start {
  my $self = shift;
  $self->{gene2_start} = shift if @_;
  return($self->{gene2_start});
}

sub gene1_end {
  my $self = shift;
  $self->{gene1_end} = shift if @_;
  return($self->{gene1_end});
}

sub gene2_end {
  my $self = shift;
  $self->{gene2_end} = shift if @_;
  return($self->{gene2_end});
}

sub exon1_num {
  my $self = shift;
  $self->{exon1_num} = shift if @_;
  return($self->{exon1_num});
}

sub exon2_num {
  my $self = shift;
  $self->{exon2_num} = shift if @_;
  return($self->{exon2_num});
}

sub exon1_id {
  my $self = shift;
  $self->{exon1_id} = shift if @_;
  return($self->{exon1_id});
}

sub exon2_id {
  my $self = shift;
  $self->{exon2_id} = shift if @_;
  return($self->{exon2_id});
}

sub transcript1_id {
  my $self = shift;
  $self->{transcript1_id} = shift if @_;
  return($self->{transcript1_id});
}

sub transcript2_id {
  my $self = shift;
  $self->{transcript2_id} = shift if @_;
  return($self->{transcript2_id});
}

sub distance1 {
  my $self = shift;
  $self->{distance1} = shift if @_;
  return($self->{distance1});
}

sub distance2 {
  my $self = shift;
  $self->{distance2} = shift if @_;
  return($self->{distance2});
}

sub format_bedpe_line {
  my ($self, $type) = @_;

  my @pe_fields = ($self->{'breakpoint'},$self->{'strand1'},$self->{'strand2'},$self->{'gene1'},$self->{'gene1_id'},$self->{'gene2'},$self->{'gene2_id'});
  if($type eq 'gene'){
    my @gene_fields = ($self->{'chr1'},$self->{'gene1_start'},$self->{'gene1_end'},$self->{'chr2'},$self->{'gene2_start'},$self->{'gene2_end'}, 'gene');
    unshift @pe_fields, @gene_fields;
  }
  elsif($type eq 'exon'){
    my @exon_fields = ($self->{'chr1'},$self->{'feature1_start'},$self->{'feature1_end'},$self->{'chr2'},$self->{'feature2_start'},$self->{'feature2_end'}, 'exon');
    unshift @pe_fields, @exon_fields;
    push @pe_fields, ($self->{'exon1_num'}, $self->{'exon1_id'},$self->{'transcript1_id'},$self->{'exon2_num'}, $self->{'exon2_id'},$self->{'transcript2_id'},$self->{'distance1'},$self->{'distance2'});
  }
  else{
    die "Only gene and exon bedpe files are currently formatted\n";
  }
  
  my $bedpe_line = join("\t", @pe_fields);
  return $bedpe_line;
}
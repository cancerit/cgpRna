package Sanger::CGP::Star::Implement;
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
use warnings FATAL => 'all';
use autodie qw(:all);
use English qw( -no_match_vars );
use Const::Fast qw(const);
use File::Which qw(which);
use File::Spec;
use File::Path qw(remove_tree make_path);
use File::Basename;
use File::Copy;
use Cwd;
use FindBin qw($Bin);
use Capture::Tiny qw(capture capture_stdout);
use PCAP::Cli;
use PCAP::Threaded;
use PCAP::Bwa::Meta;
use PCAP::Bam;
use Sanger::CGP::CgpRna;

use Data::Dumper;

const my $BAMFASTQ => q{ exclude=QCFAIL,SECONDARY,SUPPLEMENTARY T=%s S=%s O=%s O2=%s gz=1 level=1 F=%s F2=%s filename=%s};
const my $BAMFASTQ_ANALYSIS => q{ exclude=SECONDARY,SUPPLEMENTARY T=%s S=%s O=%s O2=%s gz=1 level=1 F=%s F2=%s filename=%s};
const my $FUSIONS_FILTER => q{ -i %s -s %s -n %s -o %s -p star};
const my $STAR_MAX_CORES => 16;
const my $STAR_DEFAULTS_SECTION => 'star-parameters';
const my $STAR_FUSION_SECTION => 'star-fusion-parameters';
const my $STAR => q{ %s %s --readFilesIn %s };
const my $STAR_FUSION => q{ %s --chimeric_out_sam %s --chimeric_junction %s --ref_GTF %s --min_novel_junction_support 10 --min_alt_pct_junction 10.0 --out_prefix %s };
const my $SAMTOBAM => q{ view -bS %s > %s };
const my $BAMSORT => q{ I=%s fixmate=1 inputformat=bam level=1 tmpfile=%s/tmp O=%s inputthreads=%s outputthreads=%s};
const my $HD_LINE => '@HD	VN:1.4	SO:unsorted';

sub check_input {
	my $options = shift;

  my $fusion_mode;
	
	if(exists $options->{'fusion_mode'}){
	  $fusion_mode = $options->{'fusion_mode'};
	}

	my $ref_data = $options->{'refdataloc'};
	my $species = $options->{'species'};
	my $ref_build = $options->{'referencebuild'};
	my $gene_build = $options->{'genebuild'};
	my $ref_build_loc = File::Spec->catdir($ref_data, $species, $ref_build);

	# Check the gtf and normal fusions files exist
	PCAP::Cli::file_for_reading('gtf-file', File::Spec->catfile($ref_build_loc, 'star', $gene_build, $options->{'gtffilename'}));

	if($fusion_mode){
	  PCAP::Cli::file_for_reading('normals-list',File::Spec->catfile($ref_build_loc,'cgpRna',$options->{'normalfusionslist'}));
	}
	
	my $input_meta = PCAP::Bwa::Meta::files_to_meta($options->{'tmp'}, $options->{'raw_files'}, $options->{'sample'});
	
	$options->{'meta_set'} = $input_meta;
	
	# Check the input data.
	$options->{'max_split'} = scalar @{$input_meta};
	
	my $fq = 0;
	my $fqgz = 0;
	my $bam = 0;
	
	# If the input is fastq need to check whether there is a mixture of gzipped and uncompressed fastqs.
	
	for my $input(@{$input_meta}) {
		if($input->fastq) {
			if($input->fastq =~ m/\.gz$/) {
				$fqgz = 1;
			}
			else{
				$fq = 1;
			}
		}
		else{
			$bam = 1;
		}
	}
	
	# If BAMs have been detected set flag so that they will be converted to pairs of fastq files in the prepare subroutine
	if($bam == 1){
		$options->{'bam'} = 1;
	}
	else{
		# If a mixture of gzipped and uncompressed fastqs, set flag to gzip the plain fastq/fq files in the prepare subroutine
		if($fq == $fqgz){
			$options->{'mixedfq'} = 1;
		}
	}
	
	return 1;
}

sub filter_fusions {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $sample = $options->{'sample'};
	my $star_outdir = File::Spec->catdir($options->{'tmp'}, 'star');
	my $fusions_file = File::Spec->catfile($star_outdir, "$sample.fusion_candidates.txt");
	die "The star fusion output files are missing, please run the starfusion step prior to filter.\n" unless(-e $fusions_file);

	my $normals_file = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'},'cgpRna',$options->{'normalfusionslist'});

	my $command = "$^X ";
	$command .= _which('filter_fusions.pl');
	$command .= sprintf $FUSIONS_FILTER, 	$fusions_file,
						$sample,
						$normals_file,
						$options->{'outdir'};

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	return 1;
}

sub format_rg_tags {
  my $options = shift;

  my $sample = $options->{'sample'};
  
  # Format the PU RG tag if the npg_run and lane_pos parameters have been provided
  $options->{'PU'} = $options->{'npg'}."_".$options->{'lane_pos'} if(defined $options->{'npg'} && $options->{'lane_pos'});
  
  # Check the input data
	my $input_meta = $options->{'meta_set'};
	
	# Get the RG header information to format the @RG line for the mapped BAM
	my $first_file = $input_meta->[0];
	my $rg_line;
  
  # Retrieve RG tag information for the input fastq or BAM
	if($first_file->fastq) {
		$rg_line = $first_file->rg_header(q{\t});
	}
	else {
		($rg_line, undef) = PCAP::Bam::rg_line_for_output($first_file->in, $sample);
		$rg_line = $rg_line;
	}
	
	# Prepare the new RG tags for the CGP mapped BAM
	my @rg_tags = split(/\\t/, $rg_line);
	my @comment_rg_tags = '@CO';
	my @rg_header;
	
	# Need to make sure the ID tag is the first
	if(exists $options->{'ID'}){
	  push @rg_header, 'ID:'.$options->{'ID'};
	}
	else{
	  foreach my $r(@rg_tags){
	    my @tag = split ':', $r;
	    if($tag[0] eq 'ID'){
	      push @rg_header, $r;
	      push @comment_rg_tags, 'original_'.$r;
	    }
	  }
	}

	push @rg_header, 'LB:'.$options->{'LB'} if(exists $options->{'LB'});
	# Quotes need to be around the description (DS:) tag text
	push @rg_header, 'DS:'.$options->{'DS'} if(exists $options->{'DS'});
	push @rg_header, 'PL:'.$options->{'PL'} if(exists $options->{'PL'});
	push @rg_header, 'PU:'.$options->{'PU'} if(exists $options->{'PU'});

	foreach my $r(@rg_tags){
	  unless($r eq '@RG' || $r =~ /^ID/){
	    my @tag = split ':', $r;
	    if(!exists $options->{$tag[0]}){
	      push @rg_header, $r;
	    }
	    # Once the RG tag has been formatted correctly for the CGP mapped BAM, add any pre-existing tags to a comment line to store what was in the BAM RG tags previously
	    $r = 'original_'.$r;
	    push @comment_rg_tags, $r;
	  }
	}
	
	my $comment_line = join("\t",@comment_rg_tags);
	$comment_line =~ s/"//g;
	
	# Write the old RG tags to a comment line in a file which STAR will read in using the outSAMheaderCommentFile parameter
	my $comment_file = File::Spec->catfile($options->{'tmp'}, 'star_comment_file.txt');
	open(my $ofh, '>', $comment_file) or die "Could not open file '$comment_file' $!";
  print $ofh $comment_line;
	close($ofh);
	
	$options->{'commentfile'} = $comment_file unless($first_file->fastq);

	# convert tags to shell safe tags before passing them to star in command line
	my @shell_safe_rg_header = map { _to_commandline_safe_tag_for_star($_) } @rg_header;
	$options->{'rgline'} = join(" ", @shell_safe_rg_header);

  return 1;
}

sub prepare {

	my ($index, $options) = @_;
	return 1 if(exists $options->{'index'} && $index != $options->{'index'});

	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), $index);
	
	if($options->{'bam'} || $options->{'mixedfq'}){
		my $input_meta = $options->{'meta_set'};
		my $inputdir = File::Spec->catdir($tmp, 'input');
		my $iter = 1;
		
		for my $input(@{$input_meta}) {
			next if($iter++ != $index); # skip to the relevant element in the list
			# If the file name ends fastq or fq we need to gzip it to make it consistent with other input files in the list
			if($input->fastq) {
				if($input->fastq =~ m/f(ast)?q$/) {
					my ($filename, $filepath) = fileparse($input->in);
					my $suffix = ".".$input->fastq;
					my $outfile = File::Spec->catfile($inputdir,$filename);
					my $raw_files = $options->{'raw_files'};
					
					if($input->paired_fq){

						my $infile1 = $input->in."_1".$suffix;
						my $infile2 = $input->in."_2".$suffix;
						my $outfile1 = $outfile."_1".$suffix.".gz";
						my $outfile2 = $outfile."_2".$suffix.".gz";
						
						system([0,2], "(gzip -c $infile1 > $outfile1) >& /dev/null") unless(-e $outfile1);
						system([0,2], "(gzip -c $infile2 > $outfile2) >& /dev/null") unless(-e $outfile2);
						
						for(my $i=0; $i<@$raw_files; $i++){
							@$raw_files[$i] = $outfile1 if(@$raw_files[$i] =~ $infile1);
							@$raw_files[$i] = $outfile2 if(@$raw_files[$i] =~ $infile2);
						}
					}
					else{

						my $infile = $input->in.$suffix;
						my $outfile = $outfile.$suffix.".gz";
						
						system([0,2], "(gzip -c $infile > $outfile) >& /dev/null") unless(-e $outfile);
						
						for(my $i=0; $i<@$raw_files; $i++){
							@$raw_files[$i] = $outfile if(@$raw_files[$i] =~ $infile);
						}				
					}
				}
			}
			else{
				# It's a BAM so call bamtofastq to split into gzipped paired fastq files
				my $sample = $options->{'sample'};
				my $rg = $input->{'rg'};
				
				my $fusion_mode;
	
				if(exists $options->{'fusion_mode'}){
	  		  $fusion_mode = $options->{'fusion_mode'};
				}
				
				my $command = _which('bamtofastq') || die "Unable to find 'bamtofastq' in path";
				
				if($fusion_mode){
				  # If the BAM has already been mapped and QCed and is being split for downstream analysis we just want to exclude SECONDARy and SUPPLEMENTARY reads
				  $command .= sprintf $BAMFASTQ_ANALYSIS, 	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.s"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o1"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o2"),
								File::Spec->catfile($inputdir, $sample.'.'.$rg.'_1.fastq.gz'),
								File::Spec->catfile($inputdir, $sample.'.'.$rg.'_2.fastq.gz'),
								$input->in;
				}
				else{
				  # Also exclude QCFAIL reads if the BAM is directory from iRODs or externally imported
				  $command .= sprintf $BAMFASTQ, 	File::Spec->catfile($tmp, "bamtofastq.$sample.$rg"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.s"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o1"),
								File::Spec->catfile($tmp, "bamtofastq.$sample.$rg.o2"),
								File::Spec->catfile($inputdir, $sample.'.'.$rg.'_1.fastq.gz'),
								File::Spec->catfile($inputdir, $sample.'.'.$rg.'_2.fastq.gz'),
								$input->in;
				}
																	
				PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, $index);
			
			}
		}
	}
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), $index);
	
	return 1;
}

sub process_star_params {
	my $options = shift;
	
	my $fusion_mode;
	
	if(exists $options->{'fusion_mode'}){
	  $fusion_mode = $options->{'fusion_mode'};
	}
	
	my $ini_file = $options->{'config'};
	my $cfg = new Config::IniFiles( -file => $ini_file ) or die "Could not open config file: $ini_file";
	
	die "The star default parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($STAR_DEFAULTS_SECTION));
	
	# Get the Star sample specific parameters from the options hash
	my $threads = $STAR_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $STAR_MAX_CORES);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'runThreadN', $threads);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'outFileNamePrefix', $options->{'tmp'}.'/star/');
	
	my $gtf = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'}, 'star', $options->{'genebuild'}, $options->{'gtffilename'});
	my $star_index = File::Spec->catdir($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'}, 'star');
	
	$cfg->setval($STAR_DEFAULTS_SECTION, 'sjdbGTFfile', $gtf);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'genomeDir', $star_index);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'outSAMattrRGline', $options->{'rgline'}) if(defined $options->{'rgline'});
	$cfg->setval($STAR_DEFAULTS_SECTION, 'outSAMheaderCommentFile', $options->{'commentfile'}) if(defined $options->{'commentfile'});
	$cfg->setval($STAR_DEFAULTS_SECTION, 'quantMode', 'TranscriptomeSAM') unless(defined $fusion_mode);
	$cfg->setval($STAR_DEFAULTS_SECTION, 'outSAMheaderHD', $HD_LINE);
	
	my @star_command;
	my @star_defaults = $cfg->Parameters($STAR_DEFAULTS_SECTION);
	
	for my $key(@star_defaults){
		if($cfg->val($STAR_DEFAULTS_SECTION, $key) ne ''){
			push @star_command, "--".$key." ".$cfg->val($STAR_DEFAULTS_SECTION, $key);
		}
	}
	
	if(defined $fusion_mode) {
		die "The star-fusion default parameters are missing from the config file: $ini_file\n" unless($cfg->SectionExists($STAR_FUSION_SECTION));
		
		my @star_fusion_defaults = $cfg->Parameters($STAR_FUSION_SECTION);
		for my $key(@star_fusion_defaults){
			if($cfg->val($STAR_FUSION_SECTION, $key) ne ''){
				push @star_command, "--".$key." ".$cfg->val($STAR_FUSION_SECTION, $key);
			}
		}
	}
	
	return join(" ", @star_command);
}

sub prog_version {
	my $options = shift;
	my $star_path = $options->{'starpath'};
  if(! defined $star_path || $star_path eq ''){
		$star_path = _which('STAR');
		$options->{'starpath'} = $star_path;
	}

	my $star_version;
	{
		no autodie qw(system);
		my ($stdout, $stderr, $exit) = capture{ system("$star_path --version"); };
		($star_version) = $stdout =~ /STAR_([[:digit:]\.]+.*)/m;
	}
	return $star_version;
}

sub sam_to_bam {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $sample = $options->{'sample'};
	my $outdir = $options->{'outdir'};
	my $star_outdir = File::Spec->catdir($options->{'tmp'}, 'star');
	my $command .= _which('samtools');
	$command .= sprintf $SAMTOBAM,	File::Spec->catfile($star_outdir, 'Chimeric.out.sam'),
																	File::Spec->catfile($outdir, $sample.'.star.Chimeric.out.bam');

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $command, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
	return 1;
}

sub star {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $threads = $STAR_MAX_CORES;
	$threads = $options->{'threads'} if($options->{'threads'} < $STAR_MAX_CORES);
	my $sample = $options->{'sample'};
	
	# Ensure the correct RG tags will be in the mapped BAM
	format_rg_tags($options);
	
	# Format the star command
	my $star_params = process_star_params($options);
	
	my @files1;
	my @files2;
	my $infiles1;
	my $infiles2;

	# If the input is BAM then all input fastqs will reside in the ../tmp/input directory so can search for _1 and _2 to find the files
	if($options->{'bam'}){
		my $inputdir = File::Spec->catdir($tmp, 'input');
 		opendir(my $dh, $inputdir);
		while(my $file = readdir $dh) {
			push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/_1.f/);
			push @files2, File::Spec->catfile($inputdir, $file) if($file =~ m/_2.f/);
		}
		closedir($dh);
		
		$infiles1 = join(',', sort @files1);
		$infiles2 = join(',', sort @files2);
	}
	else{
		my $raw_files = $options->{'raw_files'};
		my $input_meta = $options->{'meta_set'};
		my $first_file = $input_meta->[0];
		
		# If there are multiple input fastqs, some gzipped and some not, need to check both the input directory and raw_files array for the locations of the input files
		if($options->{'mixedfq'}){
			my $inputdir = File::Spec->catdir($tmp, 'input');
			if($first_file->paired_fq){
				for my $file(@{$raw_files}) {
					push @files1, $file if($file =~ m/_1.f.*q.gz$/);
					push @files2, $file if($file =~ m/_2.f.*q.gz$/);
				}
				opendir(my $dh, $inputdir);
				while(my $file = readdir $dh) {
					push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/_1.f.*q.gz$/);
					push @files2, File::Spec->catfile($inputdir, $file) if($file =~ m/_2.f.*q.gz$/);
				}
				closedir($dh);
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = join(',', sort @files2);
			
			}
			else{
				for my $file(@{$raw_files}) {
					push @files1, $file if($file =~ m/.f.*q.gz$/);
				}
				opendir(my $dh, $inputdir);
				while(my $file = readdir $dh) {
					push @files1, File::Spec->catfile($inputdir, $file) if($file =~ m/.f.*q.gz$/);
				}
				closedir($dh);
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = '';
			}

		}
		# Otherwise there should be just one pair of fastq files or an interleaved/single ended fastq file so only need to check the raw_files array
		else{
			if($first_file->paired_fq){
				for my $file(@{$raw_files}) {
					push @files1, $file if($file =~ m/_1.f/);
					push @files2, $file if($file =~ m/_2.f/);
				}
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = join(',', sort @files2);
				
			}
			else{
				for my $file(@{$raw_files}) {
					push @files1, $file;
				}
				
				$infiles1 = join(',', sort @files1);
				$infiles2 = '';
				
			}
		}
	}
	
	if(! defined $options->{'starpath'} || $options->{'starpath'} eq ''){
	  my $version = prog_version($options);
	}
	
	my $star_command = sprintf $STAR,	$options->{'starpath'},
					$star_params,
					$infiles1." ".$infiles2;
																
	$star_command .= "--readFilesCommand zcat" if($infiles1 =~ m/\.gz$/);
	
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $star_command, 1);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 1);
	
	

	my $stardir = File::Spec->catdir($options->{'tmp'},'star');
	my $bamsort_path = which('bamsort') || die "Unable to find 'bamsort' in path\n";
	
	my $bamsort_command1 = $bamsort_path.sprintf $BAMSORT, File::Spec->catfile($stardir, 'Aligned.out.bam'),
														$stardir,
														File::Spec->catfile($stardir, 'Aligned.sortedByCoord.out.bam'),
														$threads,
														$threads;

	my $bamsort_command2 = $bamsort_path.sprintf $BAMSORT, File::Spec->catfile($stardir, 'Aligned.toTranscriptome.out.bam'),
														$stardir,
														File::Spec->catfile($stardir, 'Aligned.toTranscriptome.sortedByCoord.out.bam'),
														$threads,
														$threads;
														
	my $fusion_mode;
	
	if(exists $options->{'fusion_mode'}){
	  $fusion_mode = $options->{'fusion_mode'};
	}
															
	my @commands;
	push @commands, $bamsort_command1;
	push @commands, $bamsort_command2 unless(defined $fusion_mode);

	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), \@commands, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
	return 1;
}

sub star_fusion {
	my $options = shift;
	
	my $tmp = $options->{'tmp'};
	return 1 if PCAP::Threaded::success_exists(File::Spec->catdir($tmp, 'progress'), 0);
	
	my $star_dir = File::Spec->catdir($tmp, 'star');
	my $chimeric_junction = File::Spec->catfile($star_dir, 'Chimeric.out.junction');
	my $chimeric_sam = File::Spec->catfile($star_dir, "Chimeric.out.sam");
	die "Some star-fusion setup files are missing, please check the output from STAR and re-run the previous stage (star) if necessary." unless(-e $chimeric_junction && -e $chimeric_sam);
	
	my $sample = $options->{'sample'};
	my $gtf = File::Spec->catfile($options->{'refdataloc'},$options->{'species'},$options->{'referencebuild'}, 'star',$options->{'genebuild'}, $options->{'gtffilename'});
	
	my $starfusionpath = $options->{'starfusionpath'};
	if(! defined $starfusionpath || $starfusionpath eq ''){
	  $starfusionpath = _which('STAR-Fusion');
	}
	
	my $commands = sprintf $STAR_FUSION,	$starfusionpath,
						$chimeric_sam,
						$chimeric_junction,
						$gtf,
						$star_dir."/$sample";
																			
	PCAP::Threaded::external_process_handler(File::Spec->catdir($tmp, 'logs'), $commands, 0);
	PCAP::Threaded::touch_success(File::Spec->catdir($tmp, 'progress'), 0);
	
	return 1;
}

sub _which {
	my $prog = shift;
	my $l_bin = $Bin;
	my $path = File::Spec->catfile($l_bin, $prog);
	$path = which($prog) unless(-e $path);
	die "Failed to find $prog in path or local bin folder ($l_bin)\n\tPATH: $ENV{PATH}\n" unless(defined $path && -e $path);
	return $path;
}

sub _to_commandline_safe_tag_for_star {
	# according to this: https://unix.stackexchange.com/a/398649
	# preserving tag values is complicated
	my $tag = shift;
	# first skip back slashes
	$tag =~ s/\\/\\\\/g;
	# then skip other needed-to-skip characters.
	# this will inevitably introduce an extra back slash into the tag if there's an '!', but for now no other way to skip '!' to prevent bash history expansion.
	my @need_to_skips = ('$', '`', '"', '!');
	foreach my $x (@need_to_skips) {
		my $pattern = '['.$x.']';
		$tag =~ s/$pattern/\\$x/g;
	}
	# double quoted it
	return '"'.$tag.'"';
}

1;

__END__

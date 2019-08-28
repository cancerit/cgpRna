import os
import sys
import re
import shutil
from string import Template
from . import run_templates_in_shell, untar, mkdir

BAMSTAT_GENOME_TEMPLATE = Template('bam_stats  -r $fai_file -i $input -o $out_dir/$sample_name.bam.bas')
BAMSTAT_TRANSCRIPTOME_TEMPLATE = Template('bam_stats  -i $trans_bam -o $out_dir/$sample_name.transcriptome.bas')
RSEQC_RRNA_TEMPLATE = Template('split_bam.py -i $input -r $ribsomal_rna_bed -o $out_dir/$sample_name.rRNA > $out_dir/$sample_name.rrna.txt')
RSEQC_GENE_COVERAGE_TEMPLATE = Template('geneBody_coverage.py -i $input -r $house_keeping_gene_bed -f png -o $out_dir/$sample_name')
RSEQC_READ_DISTRIBUTION_TEMPLATE = Template('read_distribution.py -i $input -r $reference_bed > $out_dir/$sample_name.read_dist.txt')
PROCESS_RNA_LANE_STATS_TEMPLATE = Template('process_qcstats.pl -s $sample_name -i $out_dir -o $out_dir')
COLLATE_RNA_LANE_STATS_TEMPLATE = Template('paste $out_dir/$sample_name.bam.bas $out_dir/$sample_name.insert.bas $out_dir/$sample_name.read.dist.bas $out_dir/$sample_name.rrna.bas $out_dir/$sample_name.gene.cov.bas > $out_dir/$sample_name.RNA.bas')

FAI_FILE='genome.fa.fai'
RIBSOMAL_RNA_BED='rRNA.bed'
HOUSE_KEEPING_GENE_BED='HouseKeepingGenes.bed'
REFERENCE_BED='RefSeq.bed'


def generate_stats(args):
    '''
    Top level entry point for generating stats from mapped RNA-Seq sequence files.
    '''
    # only use temp_dir when needed to extract reference files
    temp_dir = os.path.join(os.path.abspath(args.out_dir), 'cgpRna_mappingStats_temp')
    clean_temp = 0

    # guess a sample name from input file name
    # this sample_name will be used as output file name prefix.
    # NOTE: geneBody_coverage.py (from RSeQC) uses bam file name without extension as a variable name in its output r script, unless this is changed, there's no way to use other source as output file prefix.
    sample_name, _ = os.path.splitext(os.path.basename(args.input))

    # prepare the output dir
    mkdir(args.out_dir)

    # prepare the reference
    reference_data_root=os.path.abspath(args.ref)
    # If a ref bundle tar file is supplied
    # Anything else will be treated as a reference folder
    if re.search(r'\.tar\.gz$', args.ref):
        mkdir(temp_dir)
        reference_data_root=os.path.join(temp_dir, 'ref')
        untar(args.ref, reference_data_root)
    
    # gathering parameters
    params = {
        'sample_name': sample_name,
        'input': os.path.abspath(args.input),
        'out_dir': os.path.abspath(args.out_dir),
        'trans_bam': None if args.trans_bam is None else os.path.abspath(args.trans_bam),
        'fai_file': os.path.join(reference_data_root, FAI_FILE),
        'ribsomal_rna_bed': os.path.join(reference_data_root, RIBSOMAL_RNA_BED),
        'house_keeping_gene_bed': os.path.join(reference_data_root, HOUSE_KEEPING_GENE_BED),
        'reference_bed': os.path.join(reference_data_root, REFERENCE_BED)
    }

    run_templates_in_shell(
        [
            BAMSTAT_GENOME_TEMPLATE,
            RSEQC_RRNA_TEMPLATE,
            RSEQC_GENE_COVERAGE_TEMPLATE,
            RSEQC_READ_DISTRIBUTION_TEMPLATE
        ],
        params)

    if args.trans_bam is not None:
        run_templates_in_shell(
            [
                BAMSTAT_TRANSCRIPTOME_TEMPLATE, 
                PROCESS_RNA_LANE_STATS_TEMPLATE,
                COLLATE_RNA_LANE_STATS_TEMPLATE
            ],
            params)

    # clean temp dir
    if clean_temp:
        shutil.rmtree(temp_dir)

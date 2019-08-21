import os
import sys
import re
import shutil
from string import Template
from . import run_templates_in_shell, untar

BAMSTAT_GENOME_TEMPLATE=Template('bam_stats  -r $fai_file -i $input -o $out_dir/$sample_name.bam.bas')
BAMSTAT_TRANSCRIPTOME_TEMPLATE=Template('bam_stats  -i $trans_bam -o out_dir/$sample_name.transcriptome.bas')
RSEQC_RRNA_TEMPLATE=Template('split_bam.py -i $input -r $ribsomal_rna_bed -o $out_dir/$sample_name.rRNA > $out_dir/$sample_name.rrna.txt')
RSEQC_GENE_COVERAGE_TEMPLATE=Template('geneBody_coverage.py -i $input -r $house_keeping_gene_bed -f png -o $out_dir/$sample_name')
RSEQC_READ_DISTRIBUTION_TEMPLATE=Template('read_distribution.py -i $input -r $reference_bed > $out_dir/$sample_name.read_dist.txt')
PROCESS_RNA_LANE_STATS_TEMPLATE=Template('process_qcstats.pl -s $sample_name -i $out_dir -o $out_dir')
COLLATE_RNA_LANE_STATS_TEMPLATE=Template('paste $out_dir/$sample_name.bam.bas $out_dir/$sample_name.insert.bas $out_dir/$sample_name.read.dist.bas $out_dir/$sample_name.rrna.bas $out_dir/$sample_name.gene.cov.bas > $out_dir/$sample_name.RNA.bas')

FAI_FILE='genome.fa.fai'
RIBSOMAL_RNA_BED='rRNA.bed'
HOUSE_KEEPING_GENE_BED='HouseKeepingGenes.bed'
REFERENCE_BED='RefSeq.bed'

def generate_stats(args):
    '''
    Top level entry point for generating stats from mapped RNA-Seq sequence files.
    '''
    # only use temp_dir when needed to extract reference files
    temp_dir = os.path.join(args.out_dir, 'cgpRna_mappingStats_temp')
    clean_temp = 0

    # guess a sample name from input file name
    if not os.path.isfile(args.input):
        sys.exit('Error: Input: %s should be a existing file' % args.input)
    sample_name, ext = os.path.splitext(os.path.basename(args.input))
    if ext != 'bam':
        sys.exit('Error: Input: %s should be a BAM file' % args.input)

    # prepare the reference
    reference_data_root=args.ref
    # If a ref bundle tar file is supplied
    # Anything else will be treated as a reference folder
    if re.search(r'\.tar\.gz$', args.ref):
        # mkdir
        os.mkdir(temp_dir)
        reference_data_root=os.path.join(temp_dir, 'ref')
        untar(args.ref, reference_data_root)
    
    # gathering parameters
    params = {
        **vars(args),
        sample_name: sample_name,
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

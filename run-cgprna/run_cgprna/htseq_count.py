import os
import shutil
from string import Template
from . import run_templates_in_shell, untar, mkdir

BAMCLOLLATE_TEMPLATE = Template('bamcollate2 collate=1 filename=$input inputformat=bam outputformat=bam level=1 exclude=SECONDARY,SUPPLEMENTARY O=$temp_dir/tmpCollated.bam')
HTSEQ_COUNT_TEMPLATE = Template('htseq-count --format=bam --order=name --stranded="no" --type="exon" --idattr="gene_id" --mode="union" --quiet $temp_dir/tmpCollated.bam $ref | bgzip -c > $out_dir/rna_htseqcount.gz')


def count(args):
    '''
    Top level entry point for generating gene counts from mapped RNA-Seq sequence files.
    '''
    # temp_dir is for the temp bam
    temp_dir = os.path.join(os.path.abspath(args.out_dir), 'cgpRna_count_temp')

    # prepare the output dir and temp dir
    mkdir(args.out_dir)
    mkdir(temp_dir)
    
    # gathering parameters
    params = {
        'input': os.path.abspath(args.input),
        'out_dir': os.path.abspath(args.out_dir),
        'ref': os.path.abspath(args.ref),
        'temp_dir': temp_dir
    }

    run_templates_in_shell(
        [
            BAMCLOLLATE_TEMPLATE,
            HTSEQ_COUNT_TEMPLATE
        ],
        params)

    # clean temp dir
    shutil.rmtree(temp_dir)

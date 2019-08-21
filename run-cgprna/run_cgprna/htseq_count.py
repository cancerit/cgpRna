import os
import sys
import re
import tarfile
from string import Template
from . import run_shell_command

'''
# rna_bamcollate2+
bamcollate2 collate=1 filename=$IN_BAM inputformat=bam outputformat=bam level=1 exclude=SECONDARY,SUPPLEMENTARY O=$OUTDIR/tmpCollated.bam

# rna_htseqcount+
mkdir $OUTDIR/HTSeq
OUT_COUNT_GZ="$OUTDIR/HTSeq/${SAMPLE_NAME}.rna_htseqcount.gz"
htseq-count --format=bam --order=name --stranded="no" --type="exon" --idattr="gene_id" --mode="union" --quiet $OUTDIR/tmpCollated.bam ${HTSEQ_GTF} | bgzip -c > ${OUT_COUNT_GZ}

rm -r $OUTDIR/tmpCollated.bam
'''

def count(args):
    '''
    Top level entry point for generating gene counts from mapped RNA-Seq sequence files.
    '''
    print(args)

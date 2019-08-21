import os
import sys
import re
import shutil
from string import Template
from . import run_templates_in_shell, untar

STAR_MAP_TEMPLATE=Template('star_mapping.pl -s $sample_name -o $out_dir -t $threads -r $reference_data_root -sp $species -rb $ref_build -gb $gene_build -g $gene_build_gtf_name $other_options $raw_reads_string')
MARK_DUPS_TEMPLATE=Template('bammarkduplicates2 I=$out_dir/$sample_name.star.Aligned.out.bam O=$out_dir/$sample_name.bam md5=1 index=1 markthreads=$threads md5filename=$out_dir/$sample_name.bam.md5 indexfilename=$out_dir/$sample_name.bam.bai M=$out_dir/$sample_name.bam.met tmpfile=$out_dir/biormdup')
BAM_INDEX_TEMPLATE=Template('bamindex < $out_dir/$sample_name.star.AlignedtoTranscriptome.out.bam > $out_dir/$sample_name.star.AlignedtoTranscriptome.out.bam.bai')

def map_seq_files(args):
    '''
    Top level entry point for mapping RNA-Seq sequence files.
    '''
    # only use temp_dir when needed to extract reference files
    temp_dir = os.path.join(args.out_dir, 'cgpRna_map_temp')
    clean_temp = 0

    other_options = []
    if args.rg_id_tag:
        other_options.append('-lane-id %s' % args.rg_id_tag)
    if args.lb_tag:
        other_options.append('-library %s' % args.lb_tag)
    if args.ds_tag:
        other_options.append('-ds-tag %s' % args.ds_tag)
    if args.pl_tag:
        other_options.append('-machine-type %s' % args.pl_tag)
    if args.pu_tag:
        other_options.append('-npg-run %s' % args.pu_tag)

    # prepare the reference
    reference_data_root=args.ref
    # If a ref bundle tar file is supplied
    # Anything else will be treated as a reference root folder
    if re.search(r'\.tar\.gz$', args.ref):
        # mkdir
        os.mkdir(temp_dir)
        clean_temp = 1
        reference_data_root=os.path.join(temp_dir, 'ref')
        untar(args.ref, reference_data_root)
    
    # gathering parameters
    params = {
        **vars(args),
        'raw_reads_string': ' '.join(args.input),
        'reference_data_root': reference_data_root,
        'other_options': ' '.join(other_options)
    }

    run_templates_in_shell(
        [
            STAR_MAP_TEMPLATE,
            MARK_DUPS_TEMPLATE,
            BAM_INDEX_TEMPLATE
        ],
        params)
    
    # clean temp dir
    if clean_temp:
        shutil.rmtree(temp_dir)

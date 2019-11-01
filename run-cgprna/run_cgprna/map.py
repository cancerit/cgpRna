import os
import sys
import re
import shutil
import fnmatch
import copy
from string import Template
from . import run_templates_in_shell, untar, mkdir

STAR_MAP_TEMPLATE = Template('star_mapping.pl -s $sample_name -o $out_dir -t $threads -r $reference_data_root -sp $species -rb $ref_build -gb $gene_build -g $gene_build_gtf_name $other_options $raw_reads_string')
MARK_DUPS_TEMPLATE = Template('bammarkduplicates2 I=$out_dir/$sample_name.star.Aligned.out.bam O=$out_dir/$sample_name.bam md5=1 index=1 markthreads=$threads md5filename=$out_dir/$sample_name.bam.md5 indexfilename=$out_dir/$sample_name.bam.bai M=$out_dir/$sample_name.bam.met tmpfile=$out_dir/biormdup')
BAM_INDEX_TEMPLATE = Template('bamindex < $out_dir/$sample_name.star.AlignedtoTranscriptome.out.bam > $out_dir/$sample_name.star.AlignedtoTranscriptome.out.bam.bai')
RENAME_OUTPUT_TEMPLATE = Template('mv "$out_dir/${sample_name}.$file_ext" "$out_dir/${out_file_prefix}.$file_ext"')

# only because star_mapping.pl will try to find files in a particular structure
REF_RELATED_DEFAULTS = {
    'species': 'unspecified_species',
    'ref_build': 'unspecified_ref_build',
    'gene_build': 'ensembl',
    'gene_build_gtf_name': 'ensembl.gtf'
}

def map_seq_files(args):
    '''
    Top level entry point for mapping RNA-Seq sequence files.
    '''
    # args to dict to allow updates later
    args_dict = copy.deepcopy(vars(args))

    # only use temp_dir when needed to extract reference files
    temp_dir = os.path.join(os.path.abspath(args.out_dir), 'cgpRna_map_temp')
    clean_temp = 0

    other_options = []
    if args.rg_id_tag:
        other_options.append('-lane-id "%s"' % args.rg_id_tag)
    if args.lb_tag:
        other_options.append('-library "%s"' % args.lb_tag)
    if args.ds_tag:
        other_options.append('-ds-tag "%s"' % args.ds_tag)
    if args.pl_tag:
        other_options.append('-machine-type "%s"' % args.pl_tag)
    if args.pu_tag:
        other_options.append('-npg-run "%s"' % args.pu_tag)

    # prepare the output dir
    mkdir(args.out_dir)

    reference_data_root = os.path.abspath(args.ref)
    # Anything not a file will be treated as a reference root folder
    if not os.path.isfile(reference_data_root):
        if not os.path.exists(reference_data_root):
            sys.exit('Error: cound not locate directory: %s' % reference_data_root)
        if any(args_dict[ele] is None for ele in REF_RELATED_DEFAULTS.keys()):
            sys.exit(
                'Error: missing required input. When "--reference" is not a reference bundle tar file, you have to provide: %s' % ', '.join(
                    [ '--' + key.replace('_', '-') for key in REF_RELATED_DEFAULTS.keys() if args_dict[key] is None])
            )
    elif not os.path.basename(reference_data_root).endswith('.tar.gz'):
    # check if input ref file has valid file extensions
        sys.exit('Error: wrong input format. "--reference" can only be a tar.gz file or a folder.')
    else:
        # If a pre-built ref bundle tar file is supplied, prepare the reference
        mkdir(temp_dir)
        clean_temp = 1

        # use the tar file name as the ref_root_dir_name, so that in BAM header, people can tell which bundle file was used
        ref_root_dir_name = re.match(r'(.*)\.tar\.gz$', os.path.basename(reference_data_root)).group(1)
        reference_data_root=os.path.join(temp_dir, ref_root_dir_name)
        
        # set ref related values for star_mapping.pl
        for arg_name in REF_RELATED_DEFAULTS.keys():
            if args_dict[arg_name] is None:
                print('Set "%s" to default.' % arg_name)
                args_dict[arg_name] = REF_RELATED_DEFAULTS[arg_name]
            if arg_name == 'gene_build':
                print('Make sure a folder named "%s" exists in the ref bundle.' % args_dict[arg_name])
            if arg_name == 'gene_build_gtf_name':
                print('Make sure a file named: "%s" exists in the gene build folder in the bundle.' % args_dict[arg_name])

        # make the folder structure
        bundle_decompress_path = os.path.join(temp_dir, ref_root_dir_name, args_dict['species'], args_dict['ref_build'])
        # dump reference bundle
        untar(args.ref, bundle_decompress_path)

    # gathering parameters
    params = {
        **vars(args),
        'raw_reads_string': ' '.join([os.path.abspath(path) for path in args.input]),
        'reference_data_root': reference_data_root,
        'other_options': ' '.join(other_options),
        'out_dir': os.path.abspath(args.out_dir),  # overwrite the value in args with absolute path
        'species': args_dict['species'],  # overwrite the value in args
        'ref_build': args_dict['ref_build'],  # overwrite the value in args
        'gene_build': args_dict['gene_build'],  # overwrite the value in args
        'gene_build_gtf_name': args_dict['gene_build_gtf_name']  # overwrite the value in args
    }

    run_templates_in_shell(
        [
            STAR_MAP_TEMPLATE,
            MARK_DUPS_TEMPLATE,
            BAM_INDEX_TEMPLATE
        ],
        params)

    if args.out_file_prefix:
        to_rename = [
            "bam",
            "bam.bai",
            "bam.md5",
            "bam.met",
            "star.Aligned.out.bam",
            "star.AlignedtoTranscriptome.out.bam",
            "star.AlignedtoTranscriptome.out.bam.bai"
        ]
        for file_ext in to_rename:
            run_templates_in_shell(
                [
                    RENAME_OUTPUT_TEMPLATE
                ],
                {
                    'out_dir': os.path.abspath(args.out_dir),
                    'sample_name': args.sample_name,
                    'out_file_prefix': args.out_file_prefix,
                    'file_ext': file_ext
                }
            )

    # clean temp dir
    if clean_temp:
        shutil.rmtree(temp_dir)

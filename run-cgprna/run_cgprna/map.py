import os
import sys
import re
import shutil
import fnmatch
from string import Template
from . import run_templates_in_shell, untar, mkdir

STAR_MAP_TEMPLATE = Template('star_mapping.pl -s $sample_name -o $out_dir -t $threads -r $reference_data_root -sp $species -rb $ref_build -gb $gene_build -g $gene_build_gtf_name $other_options $raw_reads_string')
MARK_DUPS_TEMPLATE = Template('bammarkduplicates2 I=$out_dir/$sample_name.star.Aligned.out.bam O=$out_dir/$sample_name.bam md5=1 index=1 markthreads=$threads md5filename=$out_dir/$sample_name.bam.md5 indexfilename=$out_dir/$sample_name.bam.bai M=$out_dir/$sample_name.bam.met tmpfile=$out_dir/biormdup')
BAM_INDEX_TEMPLATE = Template('bamindex < $out_dir/$sample_name.star.AlignedtoTranscriptome.out.bam > $out_dir/$sample_name.star.AlignedtoTranscriptome.out.bam.bai')
RENAME_OUTPUT_TEMPLATE = Template('mv "$out_dir/${sample_name}.$file_ext" "$out_dir/${out_file_prefix}.$file_ext"')


def map_seq_files(args):
    '''
    Top level entry point for mapping RNA-Seq sequence files.
    '''
    # keys should be the same as what they have in command_line.py
    ref_related_args = {
        '--species': args.species,
        '--reference-build': args.ref_build,
        '--gene-build': args.gene_build,
        '--gene-build-gtf-name': args.gene_build_gtf_name
    }

    # only because star_mapping.pl will try to find files in a particular structure
    ref_related_defaults = {
        '--species': 'unspecified_species',
        '--reference-build': 'unspecified_ref_build',
        '--gene-build': 'unspecified_gene_build'
    }

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
    if not os.path.isfile(reference_data_root) and any(ele is None for ele in ref_related_args.values()):
        sys.exit(
            'Error: missing required input. When "--reference" is not a reference bundle tar file, you have to provide: %s' % ', '.join(
                [ key for key, value in ref_related_args.items() if value is None])
        )

    # check if input ref file has valid file extensions
    if not os.path.basename(reference_data_root).endswith('.tar.gz'):
        sys.exit('Error: wrong input format. "--reference" can only be a tar.gz file or a folder.')

    # If a pre-built ref bundle tar file is supplied, prepare the reference
    if re.match(r'.*\.tar\.gz$', os.path.basename(reference_data_root)):
        mkdir(temp_dir)
        clean_temp = 1

        # use the tar file name as the ref_root_dir_name, so that in BAM header, people can tell which bundle file was used
        ref_root_dir_name = re.match(r'(.*)\.tar\.gz$', os.path.basename(reference_data_root)).group(1)
        reference_data_root=os.path.join(temp_dir, ref_root_dir_name)
        
        # set ref related values for star_mapping.pl
        for key,value in ref_related_args.items():
            # GTF file name will be the same as in the bundle
            if key!= '--gene-build-gtf-name':
                if value is None:
                    ref_related_args[key] = ref_related_defaults[key]
            else:
                if value is not None:
                    print('Warning: provided "--gene-build-gtf-name" will be overwritten by the GTF file name in the reference bundle.')

        # make the folder structure
        bundle_decompress_path = os.path.join(temp_dir, ref_root_dir_name, ref_related_args['--species'], ref_related_args['--reference-build'], 'star')
        final_gtf_folder = os.path.join(bundle_decompress_path, ref_related_args['--gene-build'])
        mkdir(final_gtf_folder)

        # dump reference bundle
        untar(args.ref, bundle_decompress_path)

        # find the GTF file
        gtfs = find('*.gtf', bundle_decompress_path)
        if len(gtfs) == 1:
            ref_related_args['--gene-build-gtf-name'] = os.path.basename(gtfs[0])
            # link the file to final_gtf_folder
            os.symlink(
                gtfs[0],
                os.path.join(final_gtf_folder, ref_related_args['--gene-build-gtf-name'])
            )
        else:
            sys.exit('Error: none or too many GTF files in refence bundle. Found GTF(s): %s' % ','.join(gtfs))

    # gathering parameters
    params = {
        **vars(args),
        'raw_reads_string': ' '.join([os.path.abspath(path) for path in args.input]),
        'reference_data_root': reference_data_root,
        'other_options': ' '.join(other_options),
        'out_dir': os.path.abspath(args.out_dir),  # overwrite the value in args with absolute path
        'species': ref_related_args['--species'],  # overwrite the value in args
        'ref_build': ref_related_args['--reference-build'],  # overwrite the value in args
        'gene_build': ref_related_args['--gene-build'],  # overwrite the value in args
        'gene_build_gtf_name': ref_related_args['--gene-build-gtf-name']  # overwrite the value in args
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


def find(pattern, path):
    return (
        os.path.join(root, name)
        for root, _, files in os.walk(path) if files
        for name in files if fnmatch.fnmatch(name, pattern)
    )

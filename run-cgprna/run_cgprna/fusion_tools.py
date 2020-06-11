import os
import sys
import re
import shutil
import copy
from string import Template
from . import run_templates_in_shell, untar, mkdir
import gzip

BAM_TO_FASTQ = Template('bamtofastq exclude=SECONDARY,SUPPLEMENTARY T=$bam2fq_tmp S=$bam2fq_tmp_single_end O=$bam2fq_tmp_unmatched O2=$bam2fq_tmp_unmatched2 gz=1 level=1 F=$bam2fq_tmp_matched F2=$bam2fq_tmp_matched_2 filename=$in_bam')
# Defuse and its wapper defuse_fusion.pl cannot handle gzipped file.
BAM_TO_FASTQ_FOR_DEFUSE = Template('bamtofastq exclude=SECONDARY,SUPPLEMENTARY T=$bam2fq_tmp S=$bam2fq_tmp_single_end O=$bam2fq_tmp_unmatched O2=$bam2fq_tmp_unmatched2 F=$bam2fq_tmp_matched F2=$bam2fq_tmp_matched_2 filename=$in_bam')
TOPHAT_FUSION = Template('tophat_fusion.pl -s $sample_name -o $out_dir -t $threads -r $reference_data_root -sp $species -rb $ref_build -gb $gene_build $input')
STAR_FUSION = Template('star_fusion.pl -s $sample_name -o $out_dir -t $threads -r $reference_data_root -sp $species -rb $ref_build -gb $gene_build $input')
DEFUSE_FUSION = Template('defuse_fusion.pl -s $sample_name -o $out_dir -t $threads -r $reference_data_root -sp $species -rb $ref_build -gb $gene_build $input')
DEFUSE_FILTER = Template('defuse_filters.pl -s $sample_name -o $out_dir -i $out_dir/${sample_name}.defuse-fusion.normals.filtered.txt')

NORMAL_FUSION_ARG_NAME='-normals'

REF_RELATED_DEFAULTS = {
    'species': 'unspecified_species',
    'ref_build': 'unspecified_ref_build',
    'gene_build': 'ensembl'
}


def validate_input_seq_files(file_names):
    '''
    Perl wrappers: tophat_fusion.pl, star_fusion.pl and defuse.pl all support multiple input files but are restricted into a single type, i.e either BAM or FastQ files. As we're converting BAMs to fastqs before passing them to the wappers, we can easily handle a mixture of BAM and FastQ files, but needs to validate the files before starting to convert BAM to FastQ, so we don't waste time and resources on splitting a BAM when some of the FastQ files are not following the Perl wrapper's file name conventions.
    '''
    pairs = {}
    fq_name_pattern = re.compile(r'(.*)_([12])\.f(?:ast)?q(?:\.gz)?$')
    for a_file in file_names:
        if a_file.endswith('.bam'):
            continue
        if not os.path.exists(a_file):
            sys.exit('Error: can not find input file: %s' % a_file)
        a_file = os.path.abspath(a_file)
        match = fq_name_pattern.match(os.path.basename(a_file))
        if match.group(2):
            # if same mate number been spotted before or already got both mates
            pair_name, first_or_second_mate_in_a_pair = match.groups()
            if pairs.get(pair_name, 0) == int(first_or_second_mate_in_a_pair) or pairs.get(pair_name, 0) == 3:
                sys.exit('Error: Too many \'_%s\' mate files for prefix: %s. Possibly redundant file: %s' % (first_or_second_mate_in_a_pair, pair_name, a_file))
            pairs[pair_name] = pairs.get(pair_name, 0) + int(first_or_second_mate_in_a_pair)
        else:
            sys.exit('Error: File name does not follow expected coventions: %s' % a_file)
    for pair_name, sum_of_mate_numbers in pairs.items():
        # if both _1 and _2 of a pair of fastq files are given, value should be 3
        if sum_of_mate_numbers != 3:
            if sum_of_mate_numbers == 1:
                sys.exit('Error: Can not find second mate file of: %s' % pair_name)
            if sum_of_mate_numbers == 2:
                sys.exit('Error: Can not find first mate file of: %s' % pair_name)
            sys.exit('Error: Internal code logic error, sum of first and second mate number in a pair of fastq files is %d, which is not expected and not handled well in code.' % sum_of_mate_numbers)


# NOTE: So far only used if Defuse is given a gzipped fastq file
def gunzip(source_file, dest_file):
    print('unzipping %s ...' % os.path.basename(source_file), flush=True)
    with gzip.open(source_file, 'rb') as s_file, open(dest_file, 'wb') as d_file:
        shutil.copyfileobj(s_file, d_file, 65536)
    print('done.', flush=True)


def run_fusion_wrapper(args, temp_dir_name, fusion_templates, no_gzip=False):
    # args to dict to allow updates later
    args_dict = copy.deepcopy(vars(args))
    # only use temp_dir when needed to extract reference files
    temp_dir = os.path.join(os.path.abspath(args.out_dir), temp_dir_name)
    mkdir(temp_dir)

    # valide inputs before bamtofastq, otherwise it could be to late
    validate_input_seq_files(args.input)

    # prepare the output dir
    mkdir(args.out_dir)

    reference_data_root=os.path.abspath(args.ref)

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
        reference_data_root=os.path.join(temp_dir, 'ref')
        # set ref related args to defaults if not given
        for arg_name in REF_RELATED_DEFAULTS.keys():
            if args_dict[arg_name] is None:
                print('Set "%s" to default.' % arg_name)
                args_dict[arg_name] = REF_RELATED_DEFAULTS[arg_name]
            if arg_name == 'gene_build':
                print('Make sure a folder named "%s" exists in the ref bundle.' % args_dict[arg_name])
        untar(args.ref, os.path.join(reference_data_root, args_dict['species'], args_dict['ref_build']))

    input_fastqs = []
    fq_lane_name_count = 1

    # dumb Defuse perl wrapper doesn't like gzipped, so: 
    fq_suffix, bam_to_fq_template = ('fq.gz', BAM_TO_FASTQ) if not no_gzip else ('fq', BAM_TO_FASTQ_FOR_DEFUSE)

    for a_raw_file in args.input:
        if a_raw_file.endswith('.cram'):
            sys.exit('Error: CRAM file is not supported, please remove %s from input and retry.' % a_raw_file)
        # if input is not a bam, assume it's a fastq.
        if not a_raw_file.endswith('.bam'):
            if no_gzip and a_raw_file.endswith('.gz'):
                # dumb Defuse perl wrapper doesn't like gzipped, so:
                unzip_to = os.path.join(
                    temp_dir,
                    os.path.basename(a_raw_file)[:-3]  # remove the last 3 chars from the base name, which should be '.gz' 
                )
                gunzip(a_raw_file, unzip_to)
                input_fastqs.append(os.path.abspath(unzip_to))
            else:
                input_fastqs.append(os.path.abspath(a_raw_file))
            continue

        bam2fq_params = {
            'bam2fq_tmp': os.path.join(temp_dir, '%s.%s' % (args.sample_name, fq_lane_name_count)),
            'bam2fq_tmp_single_end': os.path.join(temp_dir, '%s.%s.s' % (args.sample_name, fq_lane_name_count)),
            'bam2fq_tmp_unmatched': os.path.join(temp_dir, '%s.%s.o1' % (args.sample_name, fq_lane_name_count)),
            'bam2fq_tmp_unmatched2': os.path.join(temp_dir, '%s.%s.o2' % (args.sample_name, fq_lane_name_count)),
            'bam2fq_tmp_matched': os.path.join(temp_dir, '%s.%s_1.%s' % (args.sample_name, fq_lane_name_count, fq_suffix)),
            'bam2fq_tmp_matched_2': os.path.join(temp_dir, '%s.%s_2.%s' % (args.sample_name, fq_lane_name_count, fq_suffix)),
            'in_bam': os.path.abspath(a_raw_file)
        }
        run_templates_in_shell([bam_to_fq_template], bam2fq_params)
        input_fastqs += [
            bam2fq_params['bam2fq_tmp_matched'],
            bam2fq_params['bam2fq_tmp_matched_2']
        ]
        fq_lane_name_count += 1

    # gathering parameters
    params = {
        'sample_name': args.sample_name,
        'input': ' '.join(input_fastqs),
        'out_dir': os.path.abspath(args.out_dir),
        'threads': args.threads,
        'reference_data_root': reference_data_root,
        'species': args_dict['species'],
        'ref_build': args_dict['ref_build'],
        'gene_build': args_dict['gene_build']
    }

    run_templates_in_shell(fusion_templates, params)

    # clean temp dir
    shutil.rmtree(temp_dir)


def tophat_fusion(args):
    '''
    Top level entry point for running tophat_fusion on RNA-Seq data.
    '''
    run_fusion_wrapper(args, 'cgpRna_tophat-fusion_temp', [TOPHAT_FUSION])


def star_fusion(args):
    '''
    Top level entry point for running star_fusion on RNA-Seq data.
    '''
    run_fusion_wrapper(args, 'cgpRna_star-fusion_temp', [STAR_FUSION])


def defuse(args):
    '''
    Top level entry point for running defuse_fusion on RNA-Seq data.
    '''
    run_fusion_wrapper(args, 'cgpRna_defuse_temp', [DEFUSE_FUSION, DEFUSE_FILTER], True)

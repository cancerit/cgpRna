"""
Handle the command line parsing and select the correct sub process.
"""

import argparse
import sys
import pkg_resources  # part of setuptools

from .map import map_seq_files
from .mapping_stats import generate_stats
from .htseq_count import count
from .bigwig import generate_bigwig

version = pkg_resources.require("run_cgprna")[0].version

def main():
    """
    Sets up the parser and handles triggereing of correct sub-command
    """
    common_parser = argparse.ArgumentParser('parent', add_help=False)
    common_parser.add_argument('-v', '--version',
                               action='version',
                               version='%(prog)s ' + version)

    parser = argparse.ArgumentParser(prog='run-cgprna', parents=[common_parser])

    subparsers = parser.add_subparsers(help='sub-command help')

    # mapping arguments
    parser_a = subparsers.add_parser(
        'map',
        parents=[common_parser],
        description='Map RNA-Seq reads to a reference genome',
        epilog='Input can be either bam or \'f(ast)?q(\.gz)?\'.')
    parser_a.add_argument(
        '-i', '--input', dest='input',
        metavar='FILE',
        nargs='+',
        help='An input raw bam file, or one of a pair of FastQ file. (optionally gzip compressed).',
        required=True)
    parser_a.add_argument(
        '-r', '--reference', dest='ref',
        metavar='TAR|PATH',
        help='A reference bundle tar file or the path to reference root directory.',
        required=True)
    parser_a.add_argument(
        '-s', '--sample-name', dest='sample_name',
        metavar='STR',
        help='Sample name, which will used to prefix output file names and SM tag in the BAM file header.',
        required=True)
    parser_a.add_argument(
        '-sp', '--species', dest='species',
        metavar='STR',
        help='Species name. e.g.: human',
        required=True)
    parser_a.add_argument(
        '-rb', '--reference-build', dest='ref_build',
        metavar='STR',
        help='Reference build name, e.g.: GRCh38.',
        required=True)
    parser_a.add_argument(
        '-gb', '--gene-build', dest='gene_build',
        metavar='STR',
        help='Gene build name, e.g.: 77.',
        required=True)
    parser_a.add_argument(
        '-gtf', '--gene-build-gtf-name', dest='gene_build_gtf_name',
        metavar='STR',
        help='File name of the gene build file, e.g.: ensemble77.gtf.',
        required=True)
    parser_a.add_argument(
        '-od', '--output_directory', dest='out_dir',
        metavar='DIR', default='.',
        help='Output directory. Default: current directory.',
        required=False)
    parser_a.add_argument(
        '-t', '--threads', dest='threads',
        metavar='INT', type=int, default=1,
        help='Number of threads to use.',
        required=False)
    parser_a.add_argument(
        '--rg-id-tag', dest='rg_id_tag',
        metavar='STR',
        help='Readgroup ID tag value in the output BAM. Default: 1 or taken from the input raw BAM file.',
        required=False)
    parser_a.add_argument(
        '--lb-tag', dest='lb_tag',
        metavar='STR',
        help='Sequencing library tag value in the output BAM header. Default: None or taken from the input raw BAM file.',
        required=False)
    parser_a.add_argument(
        '--ds-tag', dest='ds_tag',
        metavar='STR',
        help='Description tag value in the output BAM header. Default: None or taken from the input raw BAM file.',
        required=False)
    parser_a.add_argument(
        '--pl-tag', dest='pl_tag',
        metavar='STR',
        help='Platform tag value in the output BAM header. Default: None or taken from the input raw BAM file.',
        required=False)
    parser_a.add_argument(
        '--pu-tag', dest='pu_tag',
        metavar='STR',
        help='Platform unit tag value in the output BAM header. Default: None or taken from the input raw BAM file.',
        required=False)
    parser_a.set_defaults(func=map_seq_files)

    # create the parser for the "stats" command
    parser_b = subparsers.add_parser(
        'stats',
        parents=[common_parser],
        description='Generate mapping stats from a BAM file, with/without a BAM file in which reads were mapped to the transcriptome instead of genome.')
    parser_b.add_argument(
        '-i', '--input', dest='input',
        metavar='FILE',
        help='Input BAM file, in which reads are mapped to a reference genome (NOT transcriptome).',
        required=True)
    parser_b.add_argument(
        '-r', '--reference', dest='ref',
        metavar='TAR|PATH',
        help='A reference bundle tar file or the path to reference root directory.',
        required=True)
    parser_b.add_argument(
        '-tb', '--transcriptome-bam', dest='trans_bam',
        metavar='FILE',
        help='BAM file, in which reads are mapped to a reference transciptome (NOT genome).',
        required=False)
    parser_b.add_argument(
        '-od', '--output_directory', dest='out_dir',
        metavar='DIR', default='.',
        help='Output directory. Default: current directory.',
        required=False)
    parser_b.set_defaults(func=generate_stats)

    args = parser.parse_args()
    if len(sys.argv) > 1:
        args.func(args)
    else:
        sys.exit('\nERROR Arguments required\n\tPlease run: run-cgprna --help\n')

import os
from string import Template
from . import run_templates_in_shell, mkdir

BIGWIG_TEMPLATE = Template('bamToBw.pl -o $out_dir -t $threads -r $ref -b $input')


# NOTE: Require secondary input: the BAM index file
def generate_bigwig(args):
    '''
    Top level entry point for generating bigwig coverage files from mapped RNA-Seq sequence files.
    '''
    # prepare the output dir
    mkdir(args.out_dir)
    
    # gathering parameters
    params = {
        'threads': args.threads,
        'input': os.path.abspath(args.input),
        'out_dir': os.path.abspath(args.out_dir),
        'ref': os.path.abspath(args.ref),
    }

    run_templates_in_shell(
        [
            BIGWIG_TEMPLATE
        ],
        params)

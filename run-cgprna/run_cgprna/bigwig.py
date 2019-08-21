import os
import sys
import re
import tarfile
from string import Template
from . import run_shell_command

'''
# Per chromosome BigWigs
bamToBw.pl -o $OUTDIR  -r $REFERENCE -b $IN_BAM -p bamToBw

# Generates merged BigWig
bamToBw.pl -o $OUTDIR  -t 1 -r $REFERENCE -b $IN_BAM -p generateBw -i 1
'''

def generate_bigwig(args):
    '''
    Top level entry point for generating bigwig coverage files from mapped RNA-Seq sequence files.
    '''
    print(args)

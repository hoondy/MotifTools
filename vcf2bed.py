#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2018, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

###
### VCF to BED file
###

import argparse, sys
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE,SIG_DFL)

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='Takes VCF file as input and outputs in BED format')
parser.add_argument('-v','--vcf', help='VCF Input',required=True)
args = parser.parse_args()

###

with open(args.vcf, 'r') as f:
    for idx,line in enumerate(f.readlines()):
        if line.startswith("#"):
            continue
        else:

            tabs=line.strip().split("\t")
            if len(tabs)<5:
                print "Malformatted BED file. It requires 5 columns at minimum. Exiting."
                sys.exit(1)
            else:
                chr = tabs[0]
                start = str(int(tabs[1])-1)
                end = str(int(tabs[1]))
                ref = tabs[3]
                alt = tabs[4]
                print chr+"\t"+start+"\t"+end+"\t"+ref+"\t"+alt
#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2017, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

###
### BED to VCF file
###

import argparse, sys, subprocess

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='Takes BED file and outputs in VCF format')
parser.add_argument('-b','--bed', help='BED Input',required=True)
args = parser.parse_args()

###

print "##fileformat=VCFv4.2"
print "##INFO=<ID=VA,Number=.,Type=String,Description=\"Variant Annotation\">"
print "#CHROM"+"\t"+"POS"+"\t"+"ID"+"\t"+"REF"+"\t"+"ALT"+"\t"+"QUAL"+"\t"+"FILTER"+"\t"+"INFO"


subprocess.call("sort -k1,1h -k2,2n "+args.bed+" > "+args.bed+".sorted", shell=True)

with open(args.bed+".sorted", 'r') as f:
    for idx,line in enumerate(f.readlines()):
        tabs=line.strip().split("\t")
        if len(tabs)<5:
            print "Malformatted BED file. It requires 5 columns at minimum. Exiting."
            sys.exit(1)
        else:
            id = tabs[0]+":"+tabs[2]+"|"+tabs[3]+">"+tabs[4]
            extra = ""
            if len(tabs)>5:
                for i in range(5,len(tabs)):
                    extra=extra+","+tabs[i]
                extra = extra.strip(",")
            print tabs[0]+"\t"+tabs[2]+"\t"+id+"\t"+tabs[3]+"\t"+tabs[4]+"\t"+"."+"\t"+"."+"\t"+"VA="+extra

subprocess.call("rm "+args.bed+".sorted", shell=True)
#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2016, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

###
### Call Motif
###

import argparse, procMotif

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='Call TF Motif')
parser.add_argument('-s','--sample', help='Sample Name',required=True)
parser.add_argument('-t','--tf', help='TF Name',required=True)
parser.add_argument('-b','--bed', help='Peak BED File',required=True)
parser.add_argument('-r','--ref', help='REF Genome FASTA File',required=True)
args = parser.parse_args()

###

def main():

    sample = args.sample
    tf = procMotif.parseTF(args.tf)

    print "Sample:",sample
    print "TF:",tf

    m = procMotif.get_jaspar_motif(tf)

    # stop if not found
    if not m:
        print "WARNING: JASPAR motif NOT found"
        return None

    print "JASPAR motif ID:",m.matrix_id

    # getfasta
    print "Get FASTA"
    procMotif.getFasta(args.ref, args.bed, "tf_peak_"+sample+"_"+tf+".fa")
    print "DONE"

    # call TF motif
    print "Call TF motif"
    procMotif.callMotif(sample, tf, "tf_peak_"+sample+"_"+tf+".fa", args.bed, "tfbs_"+sample+"_"+tf+".bed")
    print "DONE"

    # sort & uniq
    print "sort & uniq"
    procMotif.sortUniq("tfbs_"+sample+"_"+tf+".bed","tfbs_"+sample+"_"+tf+"_uniq.bed")
    print "DONE"

main()
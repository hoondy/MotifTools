#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2016, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

import argparse, procMotif, subprocess, ConfigParser

### LOAD CONFIG ###

Config = ConfigParser.ConfigParser()
Config.read("config.ini")

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='Call TF Motif')
parser.add_argument('-s','--sample', help='Sample Name',required=True)
parser.add_argument('-t','--tf', help='TF Name',required=True)
parser.add_argument('-b','--bed', help='Peak BED File',required=True)
parser.add_argument('-r','--ref', help='REF Genome FASTA File',required=True)
args = parser.parse_args()

###

def getFasta(fastaFile, bedFile, outFile):
    subprocess.call(Config.get("app","bedtools")+" getfasta -fi "+fastaFile+" -bed "+bedFile+" -fo "+outFile, shell=True)
    return outFile

def sortUniq(inFile, outFile):
    subprocess.call("sort -k1,1 -k2,2n "+inFile+" | uniq > "+outFile, shell=True)
    return outFile

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
    getFasta(args.ref, args.bed, "tf_peak_"+sample+"_"+tf+".fa")
    print "DONE"

    # call TF motif
    print "Call TF motif"
    procMotif.callMotif(tf, "tf_peak_"+sample+"_"+tf+".fa", args.bed, float(Config.get("param","C_PVAL_THRESHOLD")), "tfbs_"+sample+"_"+tf+".bed")
    print "DONE"

    # sort & uniq
    print "sort & uniq"
    sortUniq("tfbs_"+sample+"_"+tf+".bed","tfbs_"+sample+"_"+tf+"_uniq.bed")
    print "DONE"

main()
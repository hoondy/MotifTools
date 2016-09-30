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

parser = argparse.ArgumentParser(description='Calculate TFBS D-score')
parser.add_argument('-s','--sample', help='Sample Name',required=True)
parser.add_argument('-t','--tf', help='TF Name',required=True)
parser.add_argument('-b','--bed', help='Peak BED File',required=True)
parser.add_argument('-v','--vcf', help='Variant VCF File',required=True)
parser.add_argument('-r','--ref', help='REF Genome FASTA File',required=True)
args = parser.parse_args()

###

def getPeakWithVariant(bedFile, vcfFile, outFile):
    subprocess.call(Config.get("app","bedtools")+" intersect -a "+bedFile+" -b "+vcfFile+" -wa -wb > "+outFile, shell=True)
    return outFile

def getFasta(fastaFile, bedFile, outFile):
    subprocess.call(Config.get("app","bedtools")+" getfasta -fi "+fastaFile+" -bed "+bedFile+" -fo "+outFile, shell=True)
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

    # intersect peak with variants
    print "Intersect peak with variants"
    getPeakWithVariant(args.bed, args.vcf, "tf_peak_"+sample+"_"+tf+".bed")
    print "DONE"

    # getfasta
    print "Get FASTA"
    getFasta(args.ref, "tf_peak_"+sample+"_"+tf+".bed", "tf_peak_"+sample+"_"+tf+".fa")
    print "DONE"

    # calculate D-score
    print "Calculate D-score"
    procMotif.dscoreAnalysis(tf, "tf_peak_"+sample+"_"+tf+".fa", "tf_peak_"+sample+"_"+tf+".bed", float(Config.get("param","C_PVAL_THRESHOLD")), "dscore_"+sample+"_"+tf+".bed")
    print "DONE"

main()
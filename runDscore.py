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
parser.add_argument('-b','--bed', help='Input Peak BED File',required=True)
parser.add_argument('-v','--vcf', help='Input Variant VCF File',required=True)
parser.add_argument('-r','--ref', help='REF Genome',required=True)
parser.add_argument('-a','--alt', help='ALT Genome',required=True)
args = parser.parse_args()

###

def getPeakWithVariant(bedFile, vcfFile, outFile):
    subprocess.call(Config.get("app","bedtools")+" intersect -a "+bedFile+" -b "+vcfFile+" -u > "+outFile, shell=True)
    return outFile

def getFasta(fastaFile, bedFile, outFile):
    subprocess.call(Config.get("app","bedtools")+" getfasta -fi "+fastaFile+" -bed "+bedFile+" -fo "+outFile, shell=True)
    return outFile

def main(bedFile):
    filename = bedFile.rstrip().split("/")[-1]
    filename_split = filename.rstrip().split("_")

    cell = filename_split[4]
    tf = procMotif.parseTF(filename_split[5])

    print "Processing",cell,"TF:",tf

    if procMotif.get_jaspar_motif(tf):

        # intersect peak with variants
        print "Intersect peak with variants"
        getPeakWithVariant(bedFile, args.vcf, "tf_peak_"+cell+"_"+tf+".bed")
        print "DONE"

        # getfasta
        print "Get FASTA"
        getFasta(args.ref, "tf_peak_"+cell+"_"+tf+".bed", "tf_peak_"+cell+"_"+tf+"_ref.fa")
        getFasta(args.alt, "tf_peak_"+cell+"_"+tf+".bed", "tf_peak_"+cell+"_"+tf+"_alt.fa")
        print "DONE"

        # calculate D-score
        print "Calculate D-score"
        procMotif.dscoreAnalysis(tf, "tf_peak_"+cell+"_"+tf+"_ref.fa", "tf_peak_"+cell+"_"+tf+"_alt.fa", float(Config.get("param","C_PVAL_THRESHOLD")), "dscore_"+cell+"_"+tf+".bed")
        print "DONE"

main(args.bed)
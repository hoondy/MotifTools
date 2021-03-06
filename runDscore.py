#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2016, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

###
### Disruption Score (d-score)
###

import argparse, procMotif

### LOAD ARGs ###

parser = argparse.ArgumentParser(description='Calculate TFBS/RBPBS D-score')

parser.add_argument('-n', '--name', help='Name', required=True)
parser.add_argument('-b', '--bed', help='BED File', required=True)

parser.add_argument('-m', '--motif', help='Motif File', required=True)
parser.add_argument('-f', '--format', help='Motif Format (default: jaspar) [Currently supported formats (case is ignored): AlignAce, MEME, MAST, TRANSFAC, pfm, jaspar, sites, ppm]', required=False, default="jaspar")

parser.add_argument('-v', '--vcf', help='Variant VCF File', required=True)
parser.add_argument('-r', '--ref', help='REF Genome FASTA File', required=True)

args = parser.parse_args()

###

def main():

    name = args.name

    # intersect peak with variants
    print "Intersect peak with variants"
    procMotif.getPeakWithVariant(args.bed, args.vcf, "BEDVAR_"+name+".bed")
    print "DONE"

    # getfasta
    print "Get FASTA"
    procMotif.getFasta(args.ref, "BEDVAR_"+name+".bed", "BEDVAR_"+name+".fa")
    print "DONE"

    # calculate D-score
    print "Calculate D-score"
    procMotif.dscoreAnalysis(name, args.motif, args.format, "BEDVAR_"+name+".bed", "BEDVAR_"+name+".fa", "D-SCORE_"+name+".bed")
    print "DONE"

    # sort & uniq
    print "sort & uniq"
    procMotif.sortUniq("D-SCORE_"+name+".bed","D-SCORE_"+name+"_uniq.bed")
    print "DONE"

main()
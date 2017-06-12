#!/usr/bin/python

__author__ = "Donghoon Lee"
__copyright__ = "Copyright 2016, Gerstein Lab"
__credits__ = ["Donghoon Lee"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Donghoon Lee"
__email__ = "donghoon.lee@yale.edu"

###
### Utilities
###

import argparse, procMotif, ConfigParser

### LOAD CONFIG ###

Config = ConfigParser.ConfigParser()
Config.read("config.ini")

C_PSEUDOCOUNTS = float(Config.get("param","C_PSEUDOCOUNTS"))
C_BACKGROUND_A = float(Config.get("param","C_BACKGROUND_A"))
C_BACKGROUND_C = float(Config.get("param","C_BACKGROUND_C"))
C_BACKGROUND_G = float(Config.get("param","C_BACKGROUND_G"))
C_BACKGROUND_T = float(Config.get("param","C_BACKGROUND_T"))
C_BACKGROUND = {'A':C_BACKGROUND_A,'C':C_BACKGROUND_C,'G':C_BACKGROUND_G,'T':C_BACKGROUND_T}
C_PRECISION = int(Config.get("param","C_PRECISION"))
C_PVAL_THRESHOLD = float(Config.get("param","C_PVAL_THRESHOLD"))

###

parser = argparse.ArgumentParser(description='Run Utilities')

parser.add_argument('-o','--option', help='Option',required=True)

args = parser.parse_args()

###

if args.option == "makePFM":
    procMotif.jaspar2pfm(Config.get("data","pfm_path"),Config.get("data","pfm_db_jaspar"))
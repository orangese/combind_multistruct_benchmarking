import sys
import os
import multiprocessing as mp
import wget
import ssl
import multiprocessing as mp

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def sdToMae(f):
    os.system(SCHRODINGER + "/ligprep -WAIT -adjust_itc -epik -isd " + f + " -omae " + f[:-3] + ".mae")

def parseTextMapFile(inputFile):
    os.chdir("raw_pdbs")
    with open(inputFile, 'r') as f:
        for line in f:
            line = str.strip(line)
            tempDownload = wget.download("http://files.rcsb.org/download/" + line + ".pdb", out=".")
    os.chdir("../")

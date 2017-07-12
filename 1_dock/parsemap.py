import sys
import os
import multiprocessing as mp
import wget
import ssl
import multiprocessing as mp
from schrodinger.structure import PDBWriter
from schrodinger.structure import StructureReader

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def sdToMae(f):
    os.system(SCHRODINGER + "/ligprep -WAIT -adjust_itc -epik -isd " + f + " -omae " + f[:-3] + ".mae")

def parseTextMapFile(inputFile):
    os.chdir("raw_pdbs")
    with open(inputFile, 'r') as f:
        for line in f:
            line = str.strip(line)
            #Download the PDB File from RCSB via wget
            fileName = wget.download("http://files.rcsb.org/download/" + line + ".pdb", out=".")

            #Strip out any alternate conformations
            struct = StructureReader(fileName).next()
            pw = PDBWriter(fileName, first_occ = True)
            pw.write(struct)

    os.chdir("../")

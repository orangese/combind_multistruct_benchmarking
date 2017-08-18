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

def get():
    if 'raw_pdbs' not in os.listdir('.'):
        print 'downloading pdbs...'
        input_file = [f for f in os.listdir('.') if os.path.isfile(f) and f.split('.')[-1] == 'txt'][0]
        os.system('mkdir raw_pdbs')
        with open(input_file, 'r') as f:
            for line in f:
                line = str.strip(line)
                wget.download("http://files.rcsb.org/download/" + line + ".pdb", out="./raw_pdbs")
        os.system('mv {} raw_pdbs'.format(input_file))

    print 'removing alternate conformations...'
    for f in os.listdir('raw_pdbs'):
	print str(f)
        struct = StructureReader('raw_pdbs/{}'.format(f)).next()
        pw = PDBWriter('raw_pdbs/{}.pdb'.format(f.split('.')[0].upper()), first_occ=True)
        pw.write(struct)
        if f.split('.')[0].upper() != f.split('.')[0]:
            os.system('rm raw_pdbs/{}'.format(f))

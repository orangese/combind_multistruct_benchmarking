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
    os.chdir('raw_pdbs')
    for f in [name for name in os.listdir('.') if name.split('.')[-1] == 'pdb']:
        remove_alt_conf(f)
    os.chdir('..')

def remove_alt_conf(file_name):
    struct = StructureReader(file_name).next()
    pw = PDBWriter('{}.pdb'.format(file_name.split('.')[0].upper()), first_occ=True)
    pw.write(struct)
    if file_name.split('.')[0].upper() != file_name.split('.')[0]:
        os.system('rm {}'.format(file_name))

import sys
import os
#import wget
#import ssl
from schrodinger.structure import PDBWriter, StructureReader, StructureWriter

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"

def sdToMae(f):
    os.system(SCHRODINGER + "/ligprep -WAIT -adjust_itc -epik -isd " + f + " -omae " + f[:-3] + ".mae")

def get():
#    if 'raw_pdbs' not in os.listdir('.'):
#        print 'downloading pdbs...'
#        input_file = [f for f in os.listdir('.') if os.path.isfile(f) and f.split('.')[-1] == 'txt'][0]
#        os.system('mkdir raw_pdbs')
#        with open(input_file, 'r') as f:
#            for line in f:
#                line = str.strip(line)
#                wget.download("http://files.rcsb.org/download/" + line + ".pdb", out="./raw_pdbs")
#        os.system('mv {} raw_pdbs'.format(input_file))

    print 'removing alternate conformations...'
    os.system('mkdir -p raw_maes')
    for f in os.listdir('raw_pdbs'):
        name = f.split('.')[0].upper()
        if '{}.mae'.format(name) in os.listdir('raw_maes'): 
            continue
        
        struct1 = StructureReader('raw_pdbs/{}'.format(f)).next()
        pw = PDBWriter('raw_maes/{}.pdb'.format(name, first_occ=True))
        pw.write(struct1)
        struct2 = StructureReader('raw_maes/{}.pdb'.format(name)).next()
        mw = StructureWriter('raw_maes/{}.mae'.format(name))
        mw.append(struct2)
        mw.close()
        os.system('rm raw_maes/{}.pdb'.format(name))

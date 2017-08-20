import sys
import os

def move_files():
    if 'ligands' in os.listdir('.') or 'raw_pdbs' in os.listdir('.'): 
        return
    os.system('mkdir {}/ligands'.format(s))
    os.system('mkdir {}/raw_pdbs'.format(s))
    for pdb in os.listdir('{}/original_files'.format(s)):
        os.system('cp {}/original_files/{}/{}_ligand.mol2 {}/ligands/{}_ligand.mol2'.format(s, pdb, pdb, s, pdb.upper()))
        os.system('cp {}/original_files/{}/{}_protein.pdb {}/raw_pdbs/{}.pdb'.format(s, pdb, pdb, s, pdb.upper()))

import os
import sys

from dock.chembl_props import *

def desalt(st):
    lig = None
    if len(st.molecule) == 1:
        return st
    else:
        max_atoms = []
        for mol in st.molecule:
            if len(mol.atom) > len(max_atoms):
                max_atoms = [a.index for a in mol.atom]
        return st.extract(max_atoms)

class CHEMBL:
    def __init__(self, chembl_id, smiles, ki, ki_unit, prot, val='Ki'):
        self.chembl_id = chembl_id
        self.smiles = smiles
        self.ki = ki
        self.ki_unit = ki_unit
        self.target_prot = prot
        self.val = val

        self.mw = None
        self.macrocycle = None

    def check_stereo(self):
        from schrodinger.structure import SmilesStructure
        try:
            smi_st = SmilesStructure(self.smiles).get3dStructure(True)
            self.valid_stereo = True
        except:
            smi_st = SmilesStructure(self.smiles).get3dStructure(False)
            self.valid_stereo = False

        self.st = desalt(smi_st)
        self.st.title = '{}_lig'.format(self.chembl_id)

        return self.valid_stereo

def load_chembl_proc(dir_path=None):
    ''' Read information about ligands from the chembl_info.txt file. This file is comma-delimited
    and each line in the file is formatted as so:

    chembl_id, prot, stereo, ki, smi

    Each chembl ligand is used to instantiate a CHEMBL object. In addition, this function tries to read
    molw and macrocycle data from the relevant files and link this data to the CHEMBL objects.

    Returns:
    - ligs: a dictionary from chembl ID -> CHEMBL object
    '''
    ligs = {}
    chembl_path = 'chembl/chembl_info.txt'
    if dir_path is not None:
        chembl_path = '{}/{}'.format(dir_path, chembl_path)
    if os.path.exists(chembl_path):
        with open(chembl_path) as f:
            for line in f:
                try:
                    chembl_id, prot, stereo, ki, smi = line.strip().split(',')
                    chembl_id = '{}_lig'.format(chembl_id)
                except:
                    print('chembl_info error', line)
                ligs[chembl_id] = CHEMBL(chembl_id, smi, float(ki), 'nM', prot)
                ligs[chembl_id].valid_stereo = (stereo == 'True')
    
    # Read additional ligand properties
    mw = read_molw(dir_path)
    for chembl, lig in ligs.items():
        if chembl in mw:
            lig.mw = mw[chembl]

    macrocycle = read_macrocycle(dir_path)
    for chembl, lig in ligs.items():
        if chembl in macrocycle:
            lig.macrocycle = macrocycle[chembl]
    return ligs

def load_chembl_raw(dir_path=None):
    chembl_path = 'chembl'
    if dir_path is not None:
        chembl_path = '{}/chembl'.format(dir_path)
    if not os.path.exists(chembl_path): return {}

    ligs = {}
    for c_file in os.listdir(chembl_path):
        if c_file[0] == '.' or c_file.split('.')[-1] not in ['xls', 'csv']: 
            continue
        marker = '\t' if c_file.split('.')[-1] == 'xls' else ','

        with open('{}/{}'.format(chembl_path, c_file)) as f:
            for i, line in enumerate(f):
                l_list = line.strip().split(marker)
                if not i:
                    id_ind = l_list.index('CMPD_CHEMBLID')
                    smi_ind = l_list.index('CANONICAL_SMILES')
                    type_ind = l_list.index('STANDARD_TYPE')
                    r_ind = l_list.index('RELATION')
                    val_ind = l_list.index('STANDARD_VALUE')
                    unit_ind = l_list.index('STANDARD_UNITS')
                    t_id_ind = l_list.index('PROTEIN_ACCESSION')
                    c_ind = l_list.index('CONFIDENCE_SCORE')
                    continue
                
                # Ligand criteria.
                if l_list[type_ind] not in ['Ki','IC50']: continue
                if l_list[r_ind] not in ['=']: continue
                if l_list[c_ind] not in ['9']: continue
                if l_list[unit_ind] == '': continue
                if l_list[smi_ind] == '': continue
                if l_list[unit_ind] != 'nM': continue

                cid = l_list[id_ind]
                smi = l_list[smi_ind]
                ki = float(l_list[val_ind])
                if ki <= 0: continue
                if cid not in ligs or (cid in ligs and ki < ligs[cid].ki):
                    ligs[cid] = CHEMBL(cid, smi, ki, l_list[unit_ind], 
                                       l_list[t_id_ind], l_list[type_ind])
    return ligs

def load_dude_raw(ligs, dir_path=None):
    """
    This functions very similarly to load_chembl_raw, and is meant for use with the DUD.E ligands. It
    is assumed that the ligands are packaged in .ism files in the prot/dude directory. Because the dude
    data only includes a smiles string, name, and uniprot ID (if active; decoys do not have uniprot ids),
    dummy data that does not interfere with other functions down the line is input into the CHEMBL objects.
    
    Currently, this function loads all of the avaliable actives and 100 of the decoys; this is because the
    decoy files contain thousands of ligands, but the actives may have as few as 89 (most have between 100-200).
    """
    dude_path = 'dude'
    if dir_path is not None:
        dude_path = '{}/dude'.format(dir_path)
    if not os.path.exists(dude_path): return ligs       
    
    #set dummy variables; dude data does not include these
    ki_dummy = 0.0
    unit_dummy = 'nM'
    prot_id_dummy = 12345

    for d_file in os.listdir(dude_path):
        #tagging assumes format of actives_final.ism, decoys_final.ism
        if d_file[0] == '.' or d_file.split('.')[-1] not in ['ism']:
            continue
        

        #the active dude files are in the form [smiles string] [uniprot id] [chembl name]
        #decoys take the form [smiles string] [ID]
        tag = 'active' if d_file.split('_')[0][0] == 'a' else 'decoy'              
        smi_ind = 0
        prot_id_ind = 1
        c_id_ind = 2
        if tag == 'decoy':
            c_id_ind = 1

        with open('{}/{}'.format(dude_path, d_file)) as f:
            #the dude files are in the form [smiles string] [uniprot id] [chembl name]
            for i, line in enumerate(f):
                info = line.strip().split(' ')
                smi = info[smi_ind]
                prot_id = info[prot_id_ind] if tag == 'active' else prot_id_dummy 
                cid = tag + info[c_id_ind]

                if smi == '': continue

                if cid not in ligs:
                    ligs[cid] = CHEMBL(cid, smi, ki_dummy, unit_dummy, prot_id)
                if i == 99 and tag == 'decoy': break #FOR TESTING ONLY; load 100 decoys, load all actives

    return ligs
import os
import sys

from chembl_props import *

def desalt(st):
    lig = None
    if len(st.molecule) == 1:
        return st# = smi_st
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
                ligs[chembl_id] = CHEMBL(chembl_id, smi, float(ki), 'nM', prot)#(prot, stereo, ki, smi)
                ligs[chembl_id].valid_stereo = True if stereo == 'True' else False
    
    mw = read_molw(dir_path)
    for chembl, lig in ligs.items():
        if chembl in mw:
            lig.mw = mw[chembl]

    return ligs

def load_chembl_raw(dir_path=None):
    ligs = {}

    chembl_path = 'chembl'
    if dir_path is not None:
        chembl_path = '{}/chembl'.format(dir_path)
    if not os.path.exists(chembl_path): return ligs
    for c_file in os.listdir(chembl_path):
        if c_file[0] == '.' or c_file.split('.')[-1] not in ['xls', 'csv']: 
            continue
        #print c_file
        marker = ','
        if c_file.split('.')[-1] == 'xls': marker = '\t'

        with open('{}/{}'.format(chembl_path, c_file)) as f:
            for i, line in enumerate(f):
                l_list = line.strip().split(marker)
                if i == 0:#l_list[0] == 'CMPD_CHEMBLID':# i == 0:
                    id_ind = l_list.index('CMPD_CHEMBLID')
                    smi_ind = l_list.index('CANONICAL_SMILES')
                    type_ind = l_list.index('STANDARD_TYPE')
                    r_ind = l_list.index('RELATION')
                    val_ind = l_list.index('STANDARD_VALUE')
                    unit_ind = l_list.index('STANDARD_UNITS')
                    t_id_ind = l_list.index('PROTEIN_ACCESSION')
                    c_ind = l_list.index('CONFIDENCE_SCORE')
                    #mw_ind = l_list.index('MOLWEIGHT')
                    continue
                
                if l_list[type_ind] not in ['Ki','IC50']: continue
                if l_list[r_ind] not in ['=']: continue
                #if l_list[r_ind] not in ['=','<']: continue
                if l_list[c_ind] not in ['9']: continue
                #if l_list[c_ind] not in ['8','9']: continue
                if l_list[unit_ind] == '': continue
                if l_list[smi_ind] == '': continue
                #if float(l_list[mw_ind]) > 1000: continue
 
                if l_list[unit_ind] != 'nM': continue#, l_list[id_ind]

                cid = l_list[id_ind]
                smi = l_list[smi_ind]
                ki = float(l_list[val_ind])
                if ki <= 0: continue
                if cid not in ligs or (cid in ligs and ki < ligs[cid].ki):
                    ligs[cid] = CHEMBL(cid, smi, ki, l_list[unit_ind], 
                        l_list[t_id_ind], l_list[type_ind])

    #print len(ligs)
    return ligs


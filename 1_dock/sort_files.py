import os
import sys

from schrodinger.structure import StructureReader, StructureWriter

def check_h(st):
    val = {'C':[4], 'N':[3], 'O':[2], 'S':[2,4,6], 
           'H':[1], 'P':[5], 'F':[1], 'Zn':[-2],
           'Hg':[-2], 'Cl':[1], 'K':[-1], 'Br':[1], 
           'Na':[-1], 'Ca':[-2], 'I':[1]}
    for a in st.atom:
        if a.element not in val: 
            #printa.element
            continue
        bonds = sum([b.order for b in a.bond])
        if bonds - a.formal_charge not in val[a.element]:
            print(a.element, bonds, a.index)
            return False
    return True

def split_complex(st, pdb_id):

    prot_st = st.extract([a.index for a in st.atom if a.chain != 'L'])
    prot_st.title = '{}_prot'.format(pdb_id)

    lig_path = 'structures/ligands/{}_lig.mae'.format(pdb_id)
    prot_path = 'structures/proteins/{}_prot.mae'.format(pdb_id)

    if not os.path.exists(lig_path) and len([a.index for a in st.atom if a.chain == 'L']) > 0:

        lig_st = st.extract([a.index for a in st.atom if a.chain == 'L'])
        lig_st.title = '{}_lig'.format(pdb_id)

        if not check_h(lig_st):# or not check_w(lig_st):
            print('lig error', pdb_id) 
        else:
            lig_wr = StructureWriter(lig_path)
            lig_wr.append(lig_st)
            lig_wr.close()
    
    if not os.path.exists(prot_path):
        prot_wr = StructureWriter(prot_path)
        prot_wr.append(prot_st)
        prot_wr.close()

def sort_files():#output_dir, prot):
    os.system('mkdir -p structures/proteins structures/ligands structures/complexes')

    for pdb_id in os.listdir('structures/aligned_files'):#uniprot.pdbs.items():
        opt_complex = 'structures/aligned_files/{}/{}_out.mae'.format(pdb_id, pdb_id)

        if os.path.exists(opt_complex):
            if not os.path.exists('structures/complexes/{}.mae'.format(pdb_id)):
                os.system('cp {} structures/complexes/{}.mae'.format(opt_complex, pdb_id))
            comp_st = next(StructureReader(opt_complex))
            split_complex(comp_st, pdb_id)

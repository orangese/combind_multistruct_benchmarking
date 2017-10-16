from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.structalign import StructAlign
from schrodinger.structutils.measure import get_shortest_distance
from schrodinger.structutils.analyze import AslLigandSearcher
import os

LIG_DIR = 'processed_ligands'
PROT_DIR = 'processed_proteins'

def load_proteins_and_ligands():
    proteins = {}
    ligands = {}

    for f in os.listdir(PROT_DIR):
        name = f.split('.')[0]
        proteins[name] = StructureReader('{}/{}'.format(PROT_DIR, f)).next()        
       
        lig_files = [l for l in os.listdir(LIG_DIR) if name == l.split('_')[0]] 

        if len(lig_files) == 0:
            ligands[name] = None
        else:
            ligands[name] = StructureReader('{}/{}'.format(LIG_DIR, lig_files[0])).next()
    
    return proteins, ligands

def merge_and_strip(proteins, ligands):
    merged = {}
    for struct in proteins.keys():
        all_chains = {c._getChainName() : c for c in proteins[struct].chain}
        to_delete = []
        for mol in proteins[struct].molecule:
            dist = get_shortest_distance(ligands[struct], st2=mol.extractStructure())[0]
            if dist > 10 or (len(mol.atom) == 1 and [a.element for a in mol.atom][0] == 'O'):
                to_delete.extend([a.index for a in mol.atom])
        #for c in all_chains:
        #    if c.strip() == '' or get_shortest_distance(ligands[struct], st2=all_chains[c].extractStructure())[0] > 5:
        #        to_delete.extend([a.index for a in all_chains[c].atom])
        proteins[struct].deleteAtoms(to_delete)
        merged[struct] = ligands[struct].merge(proteins[struct])
    return merged

def align():
    proteins, ligands = load_proteins_and_ligands()
    merged = merge_and_strip(proteins, ligands)
    
    os.system('rm -rf aligned_proteins aligned_ligands') 
    os.system('mkdir aligned_proteins aligned_ligands')

    structs = [s for name, s in merged.items()]

    if len(structs) > 1:
        structAlign = StructAlign()
        structAlign.align(structs[0], structs[1:])

    for struct in merged.keys():
        molecules = sorted([m for m in merged[struct].molecule], key=lambda x: len(x.residue))
        
        if ligands[struct] is not None:
            ligand = molecules[0].extractStructure()
            assert ligands[struct].isEquivalent(ligand)
            ligand._setTitle('{}_ligand'.format(struct))
            st_writer = StructureWriter('aligned_ligands/{}_ligand.mae'.format(struct))
            st_writer.append(ligand)
            st_writer.close()
            merged[struct].deleteAtoms([a.index for a in molecules[0].atom])
        st_writer = StructureWriter('aligned_proteins/{}.mae'.format(struct))
        merged[struct]._setTitle(struct)
        st_writer.append(merged[struct])
        st_writer.close()

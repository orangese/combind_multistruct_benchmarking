#!/share/PI/rondror/software/schrodinger2017-1/run
from schrodinger.structure import StructureReader
from schrodinger.structure import StructureWriter
from schrodinger.structutils.structalign import StructAlign
from schrodinger.structutils.measure import get_shortest_distance
from schrodinger.structutils.analyze import AslLigandSearcher
import os

def load_proteins_and_ligands():
    proteins = {}
    ligands = {}

    for f in os.listdir('raw_maes'):
        name = f.split('.')[0]
        proteins[name] = StructureReader('raw_maes/{}'.format(f)).next()        
       
        lig_files = [l for l in os.listdir('ligands') if name == l.split('_')[0]] 

        if len(lig_files) == 0:
            ligands[name] = None
        else:
            ligands[name] = StructureReader('ligands/{}'.format(lig_files[0])).next()
    
    return proteins, ligands

def merge_and_strip(proteins, ligands):
    merged = {}
    for struct in proteins.keys():
        all_chains = {c._getChainName() : c for c in proteins[struct].chain}
        to_delete = []
        for c in all_chains:
            if c.strip() == '' or get_shortest_distance(ligands[struct], st2=all_chains[c].extractStructure())[0] > 5:
                to_delete.extend([a.index for a in all_chains[c].atom])
        proteins[struct].deleteAtoms(to_delete)
        merged[struct] = ligands[struct].merge(proteins[struct])
    return merged

def merge_and_strip_old(proteins, ligands):

    merged = {}
    for struct in proteins.keys():
        all_chains = {c._getChainName() : c.extractStructure() for c in proteins[struct].chain if c._getChainName().strip() != ''}    
        if ligands[struct] == None:
            assert len(all_chains.keys()) == 1, 'no ligand but multiple chains for {}'.format(struct)
            merged[struct] = all_chains[all_chains.keys()[0]]
        else:
            chain_dist = [(c, get_shortest_distance(ligands[struct], st2=all_chains[c])[0]) for c in all_chains]
            chain_dist.sort(key=lambda x: x[1])
            assert chain_dist[0][1] < 3, 'no chain closer than 3 A to the ligand found in {}'.format(struct)

            if len(chain_dist) > 1:
                #assert chain_dist[1][1] > 2.5, 'chain {} is close to {}'.format(chain_dist[1][0], struct)
                #assert abs(chain_dist[0][1] - chain_dist[1][1]) > 0.5, 'chain {} and chain {} are both close to {}'.format(chain_dist[0][0], chain_dist[1][0], struct)
                if chain_dist[1][1] < 2.5 or abs(chain_dist[0][1] - chain_dist[1][1]) < 0.5: 
                    continue
            merged[struct] = ligands[struct].merge(all_chains[chain_dist[0][0]])
           
        merged[struct]._setTitle(struct)
    return merged

def align():
    proteins, ligands = load_proteins_and_ligands()
    merged = merge_and_strip(proteins, ligands)
    
    print '{} total, {} valid'.format(len(proteins.keys()), len(merged.keys()))
    os.system('rm -rf aligned_proteins aligned_ligands') 
    os.system('mkdir aligned_proteins aligned_ligands')

    structs = [s for name, s in merged.items()]

    if len(structs) > 1:
        structAlign = StructAlign()
        structAlign.align(structs[0], structs[1:])

    for struct in merged.keys():
        # ligand has one "residue"
        molecules = sorted([m for m in merged[struct].molecule], key=lambda x: len(x.residue))
        
        if ligands[struct] is not None:
            ligand = molecules[0].extractStructure()
            ligand._setTitle('{}_ligand'.format(struct))
            st_writer = StructureWriter('aligned_ligands/{}_ligand.mae'.format(struct))
            st_writer.append(ligand)
            st_writer.close()
            merged[struct].deleteAtoms([a.index for a in molecules[0].atom])
            molecules.pop(0)
        st_writer = StructureWriter('aligned_proteins/{}.mae'.format(struct))
        #for m in molecules:
            #assert len(m.residue) > 1, 'protein chain with only one residue found {}'.format(struct)
            #protein = m.extractStructure()
            #protein._setTitle(struct)
        merged[struct]._setTitle(struct)
        st_writer.append(merged[struct])
        st_writer.close()

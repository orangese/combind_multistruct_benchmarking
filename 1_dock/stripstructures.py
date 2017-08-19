#!/share/PI/rondror/software/schrodinger2017-1/run
from schrodinger.structure import StructureReader
from schrodinger.structure import StructureWriter
from schrodinger.structutils.structalign import StructAlign
from schrodinger.structutils.measure import get_shortest_distance
from schrodinger.structutils.analyze import AslLigandSearcher
import os

def load_proteins_and_ligands(protein_dir='raw_pdbs', ligand_dir='ligands'):
    
    print 'loading {} files...'.format(len(os.listdir(protein_dir)))

    proteins = {}
    ligands = {}
    for f in os.listdir(protein_dir):
        if 'pdb' not in f: continue

        name = f.split('.')[0]
        proteins[name] = StructureReader('{}/{}'.format(protein_dir, f)).next()        
        
        lig_searcher = AslLigandSearcher()
        all_ligs = lig_searcher.search(proteins[name])
        
        if ligand_dir in os.listdir('.') and len([l for l in os.listdir(ligand_dir) if name == l.split('_')[0]]) == 1:
            assert len(all_ligs) == 0, 'found both a ligand file and a ligand in the structure. delete one and try again.'
            ligands[name] = StructureReader('{}/{}'.format(ligand_dir, [l for l in os.listdir(ligand_dir) if name == l.split('_')[0]][0])).next()
            #print '+ Ligand: found in \'{}\' folder'.format(ligand_dir)
        else:
            assert len(all_ligs) in [0,1], 'multiple ligands found'
            if len(all_ligs) == 0:
                ligands[name] = None
                print '+ Ligand: no ligand found'
            if len(all_ligs) == 1:
                ligands[name] = all_ligs[0]
                #print '+ Ligand: one ligand found'
    return proteins, ligands

def merge_and_strip(proteins, ligands):

    print 'merging...'

    merged = {}
    for struct in proteins.keys():
        all_chains = {c._getChainName() : c.extractStructure() for c in proteins[struct].chain if c._getChainName().strip() != ''}    
        if ligands[struct] == None:
            assert len(all_chains.keys()) == 1, 'no ligand and multiple chains found in {}'.format(f)
            merged[struct] = all_chains[all_chains.keys()[0]]
        else:
            chain_dist = [(c, get_shortest_distance(ligands[struct], st2=all_chains[c])[0]) for c in all_chains]
            chain_dist.sort(key=lambda x: x[1])
            assert chain_dist[0][1] < 2.5, 'no chain closer than 2.5 A to the ligand found in {}'.format(f)

            if len(all_chains) > 1:
                assert abs(chain_dist[0][1] - chain_dist[1][1]) > 0.5, 'chain {} and chain {} are both close to {}'.format(chain_dist[0][0], chain_dist[1][0], struct)
                
            merged[struct] = ligands[struct].merge(all_chains[chain_dist[0][0]])
           
        merged[struct]._setTitle(struct)
    return merged

def strip(protein_dir='raw_pdbs', ligand_dir='ligands'):
    proteins, ligands = load_proteins_and_ligands(protein_dir, ligand_dir)
    merged = merge_and_strip(proteins, ligands)
    
    os.system('mkdir -p stripped ligands')
    
    print 'aligning complexes...'

    structs = [s for name, s in merged.items()]
    structAlign = StructAlign()
    structAlign.align(structs[0], structs[1:])

    for struct in merged.keys():
        # ligand has one "residue"
        molecules = sorted([m for m in merged[struct].molecule], key=lambda x: len(x.residue))
        
        if ligands[struct] is not None:
            assert len(molecules[0].residue) == 1
            assert not [r for r in molecules[0].residue][0].isStandardResidue()
            os.system('rm {}/{}*'.format(ligand_dir, struct))
            ligand = molecules.pop(0).extractStructure()
            ligand._setTitle('{}_ligand'.format(struct))
            st_writer = StructureWriter('ligands/{}_ligand.mae'.format(struct))
            st_writer.append(ligand)
            st_writer.close()
 
        st_writer = StructureWriter('stripped/{}.mae'.format(struct))
        for m in molecules:
            assert len(m.residue) > 1
            protein = m.extractStructure()
            protein._setTitle(struct)
            st_writer.append(protein)
        st_writer.close()

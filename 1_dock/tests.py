#!/share/PI/rondror/software/schrodinger2017-1/run

import sys
import os
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher

DATA = '/scratch/PI/rondror/docking_data'

PDBBIND = '/scratch/PI/rondror/docking_data/PDBbind_refined'
pdbbind_index = '{}/index/INDEX_refined_name.2016'.format(PDBBIND)
pdbbind_resolution = '{}/index/INDEX_refined_data.2016'.format(PDBBIND)

def parse_index_file():
    uniprot_to_pdb = {}
    with open(pdbbind_index) as f:
        for line in f:
            if line[0] == '#': continue
            pdb_code, year, uniprot_id = line.split()[:3]
            if uniprot_id == '------':
                uniprot_id = 'misc'
            if uniprot_id not in uniprot_to_pdb:
                uniprot_to_pdb[uniprot_id] = []
            uniprot_to_pdb[uniprot_id].append(pdb_code)
    return uniprot_to_pdb

def ready_to_move_files():
    assert 'original_files' in os.listdir('.')
    return True

def ready_to_get_pdbs():
    return True

def ready_to_extract_lig():
    return False

def get_apo():
    if 'apo.txt' not in os.listdir('.'):
        return []
    apo = []
    with open('apo.txt') as f:
        for line in f:
            apo.append(line.strip())
    return apo

def get_best_resolution():
    all_processed = [f.split('.')[0] for f in os.listdir('processed') if f.endswith('mae')]
    all_resol = {}
    with open(pdbbind_resolution) as f:
        for line in f:
            if line[0] == '#': continue
            info = line.strip().split()
            pdb = info[0].upper()
            resolution = float(info[1])
            all_resol[pdb] = resolution
    all_processed.sort(key=lambda x: all_resol[x])
    return all_processed[0]

def ready_to_align():
    assert 'ligands' in os.listdir('.'), 'ligands folder not found'
    assert 'raw_maes' in os.listdir('.'), 'raw_maes folder not found'
    
    num_pdbs = len([f for f in os.listdir('raw_pdbs') if f[:-3] == 'pdb'])
    num_maes = len([f for f in os.listdir('raw_maes') if f[:-3] == 'mae'])

    assert num_pdbs == num_maes, '# pdbs {}, # maes {}'.format(num_pdbs, num_maes)

    apo = get_apo()

    for f in os.listdir('raw_maes'):
        assert 'mae' in f, f

        name = f.split('.')[0]
        protein = StructureReader('raw_maes/{}'.format(f)).next()

        lig_searcher = AslLigandSearcher()
        all_ligs = lig_searcher.search(protein)
        assert len(all_ligs) == 0, all_ligs

        lig_files = [l for l in os.listdir('ligands') if name == l.split('_')[0]]
        if len(lig_files) == 0:
            assert name in apo, 'no ligand found for {}'.format(name)
        else:
            assert len(lig_files) == 1, lig_files
            lig = StructureReader('ligands/{}'.format(lig_files[0])).next()
            assert len([m for m in lig.molecule]) == 1, 'ligand file has too many molecules {}'.format(lig_files[0])
            res = [r for r in lig.residue]
            assert len(res) == 1, 'ligand has too many residues {}'.format(lig_files[0])
            #assert not res[0].isStandardResidue(), 'ligand should not be a standard residue {}'.format(lig_files[0])

    return True

def ready_to_process():
    assert 'aligned_ligands' in os.listdir('.')
    assert 'aligned_proteins' in os.listdir('.')
    
    #check_ligand_alignment()
    #assert len(os.listdir('aligned_ligands')) == len(os.listdir('ligands'))
    #assert len(os.listdir('aligned_proteins')) == len(os.listdir('raw_maes'))
    return True

def ready_to_check_align():
    #assert 'excluded_ligands.txt' in os.listdir('.')
    #assert 'grid_center.txt' in os.listdir('.')
    #assert 'residue_alignments.txt' in os.listdir('.')

    assert 'processed' in os.listdir('.')
    for f in os.listdir('aligned_proteins'):
        assert 'mae' in f, f
        assert f in os.listdir('processed'), f
    
    return True

def ready_to_make_grids():
    assert 'grid_center.txt' in os.listdir('.')
    #assert 'excluded_ligands.txt' in os.listdir('.')
    #assert 'residue_alignments.txt' in os.listdir('.')

    return ready_to_check_align()

def ready_to_dock():
    assert 'grids' in os.listdir('.')
    #for f in os.listdir('processed'):
    #    s = f.split('.')[0]
    #    assert s in os.listdir('grids'), s
    #    assert '{}.zip'.format(s) in os.listdir('grids/{}'.format(s)), s

    return True

def check_prerequisites(command):
    test_map = {
        'm' : ready_to_move_files,
        'r' : ready_to_get_pdbs,
        'l' : ready_to_extract_lig,
        'a' : ready_to_align,
        's' : ready_to_check_align,
        'p' : ready_to_process,
        'g' : ready_to_make_grids,
        'x' : ready_to_dock
    }
    try:
        return test_map[command]()
    except Exception as e:
        print e
        return False

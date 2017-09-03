#!/share/PI/rondror/software/schrodinger2017-1/run

import sys
import os
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher, hydrogens_present, find_overlapping_atoms
from schrodinger.structutils.measure import get_shortest_distance
from schrodinger.application.prepwizard import check_st_valences

DATA = '/scratch/PI/rondror/docking_data'

PDBBIND = '/scratch/PI/rondror/docking_data/PDBbind_refined'
pdbbind_index = '{}/index/INDEX_refined_name.2016'.format(PDBBIND)
pdbbind_info = '{}/index/INDEX_refined_data.2016'.format(PDBBIND)

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

def get_all_info():
    all_info = {}
    with open(pdbbind_info) as f:
        for line in f:
            if line[0] == '#': continue
            info = line.strip().split()
            pdb = info[0].upper()
            resolution = float(info[1])
            year = int(info[2])
            all_info[pdb] = (resolution, year)
    return all_info

def get_best_resolution():
    all_processed = sorted([f.split('.')[0] for f in os.listdir('aligned_proteins') if f.endswith('mae')])
    all_info = get_all_info()
    all_processed.sort(key=lambda x: all_info[x][0])
    for i in all_processed:
        if all_info[i] == all_resol[all_processed[0]]:
            print i, all_info[i]
    return all_processed[0]

def get_first():
    all_structs = sorted([f.upper() for f in os.listdir('original_files')])
    all_info = get_all_info()
    all_structs.sort(key=lambda x: all_info[x][1])
    for num, i in enumerate(all_structs):
        if all_info[i][1] == all_info[all_structs[0]][1]:
            print i, all_info[i]
        #assert num < 1, 'multiple structures with same year'
    return all_structs[0]

def check_metals(data_dir):
    for dataset in os.listdir(data_dir):
        print dataset
        os.chdir('{}/{}'.format(data_dir, dataset))
        count = 0
        for pdb in [s.split('.')[0] for s in os.listdir('raw_maes')]:
            lig_st = StructureReader('ligands/{}_ligand.mol2'.format(pdb)).next()
            prot_st = StructureReader('raw_maes/{}.mae'.format(pdb)).next()
            for c in prot_st.chain:
                for mol in c.extractStructure().molecule:#prot_st.molecule:
                    dist = get_shortest_distance(lig_st, st2=mol.extractStructure())[0]
                    if dist < 5 and len(mol.atom) == 1 and [a.element for a in mol.atom][0] != 'O':
                        print pdb, dist, mol.number, [a.element for a in mol.atom][0]
                        count += 1
                        break
        print '{} has {} complexes. {} have a metal present'.format(dataset, len(os.listdir('raw_maes')), count)

def check_duplicates(data_dir):
    for dataset in os.listdir(data_dir):
        print dataset
        ligands = os.listdir('{}/{}/processed_ligands'.format(data_dir, dataset))
        duplicates = {}
        all_resol = get_all_info()
        for i in range(len(ligands)):
            ri = all_resol[ligands[i].split('_')[0]][0]
            for (lig2, r2) in duplicates:
                l1 = StructureReader('{}/{}/processed_ligands/{}'.format(data_dir, dataset, ligands[i])).next()
                l2 = StructureReader('{}/{}/processed_ligands/{}'.format(data_dir, dataset, lig2)).next()
                if l1.isEquivalent(l2, False):
                    duplicates[(lig2, r2)].append((ligands[i], ri))
                    break
            else:
                duplicates[(ligands[i], ri)] = []
        
        print '{} ligands, {} unique ligands'.format(len(ligands), len(duplicates))
        
        grid = os.listdir('{}/{}/grids'.format(data_dir, dataset))[0]
        #print grid
        #os.mkdir('{}/{}/unique_ligands'.format(data_dir, dataset))
        for d in duplicates:
            identical = sorted([d] + duplicates[d], key=lambda (s,r): -1 if s == '{}_ligand.mae'.format(grid) else r)
            if grid in [i[0].split('_')[0] for i in identical]: print identical
            #os.system('cp {}/{}/processed_ligands/{} {}/{}/unique_ligands'.format(data_dir, dataset, identical[0][0], data_dir, dataset))
            #if len(identical) > 1:
            #    print identical
        #break
        #print '{} ligands, {} unique ligands'.format(len(os.listdir('{}/{}/processed_ligands'.format(data_dir, dataset))), 
        #                                             len(os.listdir('{}/{}/unique_ligands'.format(data_dir, dataset))))

def test_pdbbind(data_dir):
    all_pdbbind_refined = parse_index_file()

    for dataset in os.listdir(data_dir):
        os.chdir(data_dir)
        output_file = '{}/{}/status.txt'.format(data_dir, dataset)
        try:
            expected_structures = set([s.upper() for s in all_pdbbind_refined[dataset]])
            original_files = set([s.upper() for s in os.listdir('{}/original_files'.format(dataset))])
            assert all_files_present(expected_structures, original_files, output_file, 'index file', 'original_files')
        
            raw_pdbs = set([s.split('.')[0] for s in os.listdir('{}/raw_pdbs'.format(dataset))])
            assert all_files_present(original_files, raw_pdbs, output_file, 'original_files', 'raw_pdbs')

            ligands = set([s.split('_')[0] for s in os.listdir('{}/ligands'.format(dataset))])
            assert all_files_present(raw_pdbs, ligands, output_file, 'raw_pdbs', 'ligands')

            raw_maes = set([s.split('.')[0] for s in os.listdir('{}/raw_maes'.format(dataset))])
            assert all_files_present(raw_pdbs, raw_maes, output_file, 'raw_pdbs', 'raw_maes')
 
            aligned_ligands = set([s.split('_')[0] for s in os.listdir('{}/aligned_ligands'.format(dataset))])
            assert all_files_present(raw_maes, aligned_ligands, output_file, 'raw_maes', 'aligned_ligands')

            aligned_proteins = set([s.split('.')[0] for s in os.listdir('{}/aligned_proteins'.format(dataset))])
            assert all_files_present(raw_maes, aligned_proteins, output_file, 'raw_maes', 'aligned_proteins')

            processed_ligands = set([s.split('_')[0] for s in os.listdir('{}/processed_ligands'.format(dataset))])
            assert all_files_present(aligned_ligands, processed_ligands, output_file, 'aligned_ligands', 'processed_ligands')

            #processed_proteins = set([s.split('.')[0] for s in os.listdir('{}/processed_proteins'.format(dataset))])
            #assert all_files_present(aligned_proteins, processed_proteins, output_file, 'aligned_proteins', 'processed_proteins')
            
            os.chdir('{}/processed_ligands'.format(dataset))
            assert all_files_processed(output_file)
            os.chdir('../processed_proteins')
            assert all_files_processed(output_file)
            
            print '{} ready to go'.format(dataset)
        except AssertionError as e:
            print e
            print '{} not ready'.format(dataset)


def all_files_processed(output_file):
    passed = True
    with open(output_file, 'a') as f:
        for struct_file in sorted(os.listdir('.')):
            st = StructureReader(struct_file).next()
            if len([i for i in StructureReader(struct_file)]) > 1:
                f.write('{}: multiple structs in file\n'.format(struct_file))
                passed = False
            if st._getTitle() != struct_file.split('.')[0]:
                f.write('{}: struct name {}\n'.format(struct_file, st._getTitle()))
                passed = False
            if not hydrogens_present(st):
                f.write('{}: missing hydrogens\n'.format(struct_file))
                passed = False
            if len(find_overlapping_atoms(st, False, False)) > 0:
                f.write('{}: overlapping atoms\n'.format(struct_file))
                passed = False
            try:
                check_st_valences(st)
            except RuntimeError:
                f.write('{}: bad valences\n'.format(struct_file))
            for r in st.residue:
                if r.hasMissingAtoms():
                    f.write('{}: missing atoms\n'.format(struct_file))
                    passed = False
                    break
    return passed


def all_files_present(d1, d2, output_file, name1='', name2=''):
    passed = True
    with open(output_file, 'a') as f:
        f.write('Comparing {} and {}.\n'.format(name1, name2))
        f.write('Expected {} structures, found {}.\n'.format(len(d1), len(d2)))
        for s in d1:
            if s not in d2:
                f.write('Missing file: {}\n'.format(s))
                passed = False
        for s in d2:
            if s not in d1:
                f.write('Unexpected file: {}\n'.format(s))
                passed = False
    return passed

def ready_to_move_files():
    assert 'original_files' in os.listdir('.')
    return True

def ready_to_get_pdbs():
    return True

def ready_to_extract_lig():
    return os.path.exists('raw_maes')

def get_apo():
    if 'apo.txt' not in os.listdir('.'):
        return []
    apo = []
    with open('apo.txt') as f:
        for line in f:
            apo.append(line.strip())
    return apo

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
    return True
    assert 'processed' in os.listdir('.'), 'procesed dir not found'
    for f in os.listdir('aligned_proteins'):
        assert 'mae' in f, f
        assert f in os.listdir('processed'), f
        print f   
        st = StructureReader('processed/{}'.format(f))
        assert len([i for i in st]) == 1, '{} multiple structs in file'.format(f)
        st = st.next()
        assert st._getTitle() == f.split('.')[0], 'file {} struct name {}'.format(f, st._getTitle())
 
    return True

def ready_to_make_grids():
    assert 'grid_center.txt' in os.listdir('.'), 'grid center not found'
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

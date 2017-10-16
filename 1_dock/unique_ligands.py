import os
from schrodinger.structure import StructureReader, StructureWriter
from tests import get_all_info

def check_valence(st, st_name):
    if st._getTitle() != st_name and st._getTitle() != '{}_ligand'.format(st_name):
        print '-', st_name, st._getTitle(), 'wrong title'
    valence = {'H':set([1]), 'C':set([4]), 'N':set([3]), 'O':set([2]), 'S':set([2,4,6]), 'P':set([5])}
    for a in st.atom:
        if a.element in valence and sum([b.order for b in a.bond]) - a.formal_charge not in valence[a.element]:
            print '-', st._getTitle(), a.element, a.bond_total, a.formal_charge, a.index
            return False
    return True

def manage_ligands():
    #if process_ligands():
    filter_duplicates()


def process_proteins():
    os.system('mkdir -p broken_proteins')

    done = True
    for p in sorted(os.listdir('raw_maes')):
        st_name = p.split('.')[0].upper()
        f_name = '{}.mae'.format(st_name)
        if f_name in os.listdir('processed_proteins'):
            if not check_valence(StructureReader('processed_proteins/{}'.format(f_name)).next(), st_name):
                os.system('mv processed_proteins/{} broken_proteins'.format(f_name, f_name))
                done = False
        elif f_name in os.listdir('broken_proteins'):
            if check_valence(StructureReader('broken_proteins/{}'.format(f_name)).next(), st_name):
                os.system('mv broken_proteins/{} processed_proteins'.format(f_name))
        else:
            print '-', p, 'missing'

    if 'broken_proteins' in os.listdir('.') and len(os.listdir('broken_proteins')) == 0:
        os.system('rm -r broken_proteins')

    return done 

def process_ligands():
    os.system('mkdir -p processed_ligands broken_ligands')

    done = True
    for l in sorted(os.listdir('ligands')):
        st_name = l.split('_')[0].upper()
        f_name = '{}_ligand.mae'.format(st_name)
        if f_name in os.listdir('processed_ligands'):
            if not check_valence(StructureReader('processed_ligands/{}'.format(f_name)).next(), st_name):
                os.system('mv processed_ligands/{} broken_ligands'.format(f_name, f_name))
                done = False
        elif f_name in os.listdir('broken_ligands'):
            if check_valence(StructureReader('broken_ligands/{}'.format(f_name)).next(), st_name):
                os.system('mv broken_ligands/{} processed_ligands'.format(f_name))
        else:
            st = StructureReader('ligands/{}'.format(l)).next()
            st._setTitle('{}_ligand'.format(st_name))
            
            if not check_valence(st, st_name):
                st_writer = StructureWriter('broken_ligands/{}'.format(f_name))
                done = False
            else:
                st_writer = StructureWriter('processed_ligands/{}'.format(f_name))
            st_writer.append(st)
            st_writer.close()

    if 'broken_ligands' in os.listdir('.') and len(os.listdir('broken_ligands')) == 0:
        os.system('rm -r broken_ligands')

    return done 

def filter_duplicates():
    os.system('mkdir -p unique_ligands')

    duplicates = {}
    all_resol = get_all_info()
    for l1 in os.listdir('processed_ligands'):
        r1 = all_resol[l1.split('_')[0]][0]
        for (l2, r2) in duplicates:
            st_l1 = StructureReader('processed_ligands/{}'.format(l1)).next()
            st_l2 = StructureReader('processed_ligands/{}'.format(l2)).next()
            if st_l1.isEquivalent(st_l2, False):
                duplicates[(l2, r2)].append((l1, r1))
                break
        else:
            duplicates[(l1, r1)] = [(l1, r1)]
    print '{} ligands, {} unique ligands'.format(len(os.listdir('processed_ligands')), len(duplicates))
    for (l, r), dups in duplicates.items():
        dups.sort(key=lambda x: x[0]) # break ties alphabetically (roughly chronologically)
        dups.sort(key=lambda x: x[1]) # sort by resolution
    #    if dups[0][0] not in os.listdir('unique_ligands'):
        if len(dups) > 1:
            print dups
            for li, ri in dups[1:]:
                assert li not in os.listdir('unique_ligands'), li
            #os.system('cp processed_ligands/{} unique_ligands'.format(dups[0][0]))

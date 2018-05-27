import os
import sys

from schrodinger.structure import StructureReader, StructureWriter, _StructureProperty
from schrodinger.structutils.rmsd import *#renumber_conformer, calculate_in_place_rmsd
from schrodinger.structutils.analyze import evaluate_smarts

lig1 = sys.argv[1]
lig2 = sys.argv[2]
struct = sys.argv[3]
glide_dir = sys.argv[4]
out_dir = sys.argv[5]
mcss_dir = '../../ligands/mcss/{}'.format(out_dir)

def neutral(st,aggressive=True):
    for a in st.atom:
        a.formal_charge = 0#_setAtomFormalCharge(0)
    if aggressive:
        for b in st.bond:
            b.order = 1#_setOrder(1)
    st.retype()
    return st

def load_mcss(lig1, lig2):
    with open('{}/{}-{}.txt'.format(mcss_dir, lig1, lig2)) as f:
        l1,s1,m1,smarts1 = f.readline().strip().split(',')
        l2,s2,m2,smarts2 = f.readline().strip().split(',')
    assert lig1 == l1
    assert lig2 == l2
    assert m1 == m2
    return smarts1, smarts2, int(s1), int(s2), int(m1)

def find_mcss_matches(st1, smarts1, st2, smarts2):
    try: mcss1 = evaluate_smarts(st1, smarts1, unique_sets=True)
    except: mcss1 = evaluate_smarts(st1, smarts2, unique_sets=True)
    try: mcss2 = evaluate_smarts(st2, smarts2, unique_sets=True)
    except: mcss2 = evaluate_smarts(st2, smarts1, unique_sets=True)
    mcss1_st = [neutral(st1.extract(x)) for x in mcss1]
    mcss2_st = [neutral(st2.extract(x)) for x in mcss2]
    #print smarts1, smarts2
    #print len(mcss1), len(mcss2)
    valid_pairs = []
    for i,m1 in enumerate(mcss1_st):
        for j,m2 in enumerate(mcss2_st):
            #print sorted([a.order for a in m1.bond])
            #print sorted([a.order for a in m2.bond])
            if m1.isEquivalent(m2, False):
                valid_pairs.append((i,j))
    return mcss1, mcss2, valid_pairs

def proc_mcss(lig1, lig2, struct, glide_dir, out_dir):

    f1 = '{}-to-{}'.format(lig1, struct)
    f2 = '{}-to-{}'.format(lig2, struct)
     
    smarts1, smarts2, s1, s2, m = load_mcss(lig1, lig2)

    st1 = StructureReader('../../ligands/prepared_ligands/{}/{}.mae'.format(lig1,lig1)).next()
    st2 = StructureReader('../../ligands/prepared_ligands/{}/{}.mae'.format(lig2,lig2)).next()
    mcss1, mcss2, valid_mcss_pairs = find_mcss_matches(st1, smarts1, st2, smarts2)
   
    if len(valid_mcss_pairs) == 0:
        print 'no mcss found', lig1, smarts1, lig2, smarts2
        with open('{}-{}-to-{}.csv'.format(lig1,lig2,struct), 'w') as f:
            f.write('{},{},0\n'.format(s1, s2))
        return
    #print len(valid_mcss_pairs), 'mcss found'
    #return
    pv1 = [(x[0]-1,x[1]) for x in enumerate(StructureReader('../../{}/{}/{}_pv.maegz'.format(glide_dir, f1, f1)))][1:]
    pv2 = [(x[0]-1,x[1]) for x in enumerate(StructureReader('../../{}/{}/{}_pv.maegz'.format(glide_dir, f2, f2)))][1:]
    
    with open('{}-{}-to-{}.csv'.format(lig1,lig2,struct), 'w') as f:
        f.write('{},{},{}\n'.format(s1, s2, m))
        for i, p1 in pv1:
            for j, p2 in pv2:
                if i > 105 or j > 105: 
                    continue

                mcss1, mcss2, valid_mcss_pairs = find_mcss_matches(p1, smarts1, p2, smarts2)
                if len(valid_mcss_pairs) == 0:
                    print 'no mcss found for poses', lig1, lig2, st
                    return

                rmsd = 10000
                for m1, m2 in valid_mcss_pairs:
                    nrmsd = calculate_in_place_rmsd(p1, mcss1[m1], p2, mcss2[m2], use_symmetry=True)
                    rmsd = min(rmsd, nrmsd)

                f.write('{},{},{}\n'.format(i, j, rmsd))

proc_mcss(lig1, lig2, struct, glide_dir, out_dir)



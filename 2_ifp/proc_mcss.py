import os
import sys

from schrodinger.structure import StructureReader, StructureWriter, _StructureProperty
from schrodinger.structutils.rmsd import *#renumber_conformer, calculate_in_place_rmsd
from schrodinger.structutils.analyze import evaluate_smarts

# 1: load in reference ligand with mcss info
# 2: iterate through pose viewers

lig1 = sys.argv[1]
lig2 = sys.argv[2]
struct = sys.argv[3]
glide_dir = sys.argv[4]
out_dir = sys.argv[5]

class RefLig:
    def __init__(self, st):
        prop = _StructureProperty(st)
        self.name = st._getTitle()
        self.st = st
        self.mcss_smarts = prop['s_canvas_MCS_SMARTS']
        #self.mcss = evaluate_smarts(st, self.mcss_smarts, unique_sets=True)#, verbose=True, unique_sets=True)
        #self.mcss = [x for x in self.mcss if len(x) > 6] # ignore tiny mcss
        #self.mcss_st = [neutral(st.extract(x)) for x in self.mcss]
        self.size = len([a for a in st.atom if a.element != 'H'])

def neutral(st):
    for a in st.atom:
        a._setAtomFormalCharge(0)
    #for b in st.bond:
    #    b._setOrder(1)
    return st

def find_mcss_matches(st1, smarts1, st2, smarts2):
    mcss1 = evaluate_smarts(st1, smarts1, unique_sets=True)
    mcss2 = evaluate_smarts(st2, smarts2, unique_sets=True)
    mcss1_st = [neutral(st1.extract(x)) for x in mcss1 if len(x) > 6]
    mcss2_st = [neutral(st2.extract(x)) for x in mcss2 if len(x) > 6]
    valid_pairs = []
    for i,m1 in enumerate(mcss1_st):
        for j,m2 in enumerate(mcss2_st):
            if m1.isEquivalent(m2, False):
                valid_pairs.append((i,j))
    return mcss1, mcss2, valid_pairs

def proc_mcss(lig1, lig2, struct, glide_dir, out_dir):

    os.chdir('../../..')

    name = '{}-{}'.format(lig1, lig2)

    f1 = '{}-to-{}'.format(lig1, struct)
    f2 = '{}-to-{}'.format(lig2, struct)

    ref_st = {}
    try:
        st_out = StructureReader('ligands/mcss/{}.mae'.format(name))
        for i, st in enumerate(st_out):
            ref_st[st._getTitle()] = RefLig(st)
    except:
        print 'invalid mcss', name
        os.system('rm ligands/mcss/{}.mae'.format(name))
        return

    assert lig1 in ref_st and lig2 in ref_st and len(ref_st.keys()) == 2, ref_st.keys()
     
    smarts1, smarts2 = ref_st[lig1].mcss_smarts, ref_st[lig2].mcss_smarts

    mcss1, mcss2, valid_mcss_pairs = find_mcss_matches(ref_st[lig1].st, smarts1, ref_st[lig2].st, smarts2)
   
    if len(valid_mcss_pairs) == 0:
        with open('mcss/{}/{}/{}.csv'.format(out_dir, struct, name), 'w') as f:
            f.write('{},{},None'.format(ref_st[lig1].size, ref_st[lig2].size))
        return

    pv1 = [(x[0]-1,x[1]) for x in enumerate(StructureReader('docking/{}/{}/{}_pv.maegz'.format(glide_dir, f1, f1)))][1:]
    pv2 = [(x[0]-1,x[1]) for x in enumerate(StructureReader('docking/{}/{}/{}_pv.maegz'.format(glide_dir, f2, f2)))][1:]
    
    with open('mcss/{}/{}/{}.csv'.format(out_dir, struct, name), 'w') as f:
        f.write('{},{},{}\n'.format(ref_st[lig1].size, ref_st[lig2].size, len(mcss1[0])))
        for i, p1 in pv1:
            for j, p2 in pv2:
                if i > 105 or j > 105: 
                    continue
                #try:            
                #    renumber_conformer(ref_st[lig1].st, p1)#, use_symmetry=True)
                #    renumber_conformer(ref_st[lig2].st, p2)#, use_symmetry=True)
                #except:
                #    print i,j,len(ref_st[lig1].st.atom), len(p1.atom)                
                #    print i,j,len(ref_st[lig2].st.atom), len(p2.atom)                
                #    exit()

                mcss1, mcss2, valid_mcss_pairs = find_mcss_matches(p1, smarts1, p2, smarts2)

                #p1, p2 = neutral(p1), neutral(p2)

                rmsd = 10000
                for m1, m2 in valid_mcss_pairs:
                    nrmsd = calculate_in_place_rmsd(p1, mcss1[m1], p2, mcss2[m2], use_symmetry=True)
                    rmsd = min(rmsd, nrmsd)

                if rmsd == 10000:
                    print smarts1, smarts2, len(mcss1), len(mcss2)
                    print ref_st[lig1].st.isEquivalent(p1), ref_st[lig2].st.isEquivalent(p2)

                f.write('{},{},{}\n'.format(i, j, rmsd))#print i, j, rmsd

proc_mcss(lig1, lig2, struct, glide_dir, out_dir)



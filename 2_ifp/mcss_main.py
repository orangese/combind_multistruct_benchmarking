import os
import sys

from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.rmsd import *#renumber_conformer, calculate_in_place_rmsd
from schrodinger.structutils.analyze import evaluate_smarts

def st_reduce(st):
    heavy = st.extract([a.index for a in st.atom if a.element != 'H'])
    for b in heavy.bond:
        b.order = 1 
    for a in heavy.atom:
        a.formal_charge = 0 
    return heavy

def load_mcss(lig1, lig2, mcss_dir):
    try:
        with open('{}/{}-{}.txt'.format(mcss_dir, lig1, lig2)) as f:
            l1,s1,m1,smarts1 = f.readline().strip().split(',')
            l2,s2,m2,smarts2 = f.readline().strip().split(',')
    except:
        print 'mcss file error'
        os.system('rm -rf {}/{}-{}.txt'.format(mcss_dir, lig1, lig2))
        exit()
    assert lig1 == l1
    assert lig2 == l2
    assert m1 == m2
    return smarts1, smarts2, int(s1), int(s2), int(m1)

def find_mcss_matches(st1, smarts1, st2, smarts2):
    try: mcss1 = evaluate_smarts(st1, smarts1, unique_sets=True)
    except: mcss1 = evaluate_smarts(st1, smarts2, unique_sets=True)
    try: mcss2 = evaluate_smarts(st2, smarts2, unique_sets=True)
    except: mcss2 = evaluate_smarts(st2, smarts1, unique_sets=True)
    mcss1_st = [st1.extract(x) for x in mcss1]
    mcss2_st = [st2.extract(x) for x in mcss2]
    valid_pairs = []
    for i,m1 in enumerate(mcss1_st):
        for j,m2 in enumerate(mcss2_st):
            if m1.isEquivalent(m2, False):
                valid_pairs.append((i,j))
    return mcss1, mcss2, valid_pairs

scaff_path = '../../ligands/prepared_ligands/{}/{}_neutral.mae'
pv_path = '../../{}/{}-to-{}/{}-to-{}_pv.maegz'
outf = '{}-{}-to-{}.csv'

class MCSS:
    def __init__(self, l1, l2, st, gdir, mdir):
        self.l1 = l1
        self.l2 = l2
        self.st = st
        self.gdir = gdir
        self.mdir = mdir

        self.sm1, self.sm2, self.s1, self.s2, self.s3 = load_mcss(l1, l2, mdir)

        self.ref1 = StructureReader(scaff_path.format(l1,l1)).next()
        self.ref2 = StructureReader(scaff_path.format(l2,l2)).next()

        self.pv1 = list(StructureReader(pv_path.format(gdir, l1, st, l1, st)))[1:]
        self.pv2 = list(StructureReader(pv_path.format(gdir, l2, st, l2, st)))[1:]

        self.proc_mcss()

    def proc_mcss(self):

        with open(outf.format(self.l1,self.l2,self.st), 'w') as f:
            f.write('{},{},{}\n'.format(self.s1, self.s2, self.s3))
            if self.s3 == 0: return

            pv1 = [st_reduce(st) for st in self.pv1]
            pv2 = [st_reduce(st) for st in self.pv2]

            for st in pv1: assert st.isEquivalent(self.ref1), 'equivalence error {}'.format(self.l1)
            for st in pv2: assert st.isEquivalent(self.ref2), 'equivalence error {}'.format(self.l2)

            for i, p1 in enumerate(pv1):
                for j, p2 in enumerate(pv2):
                    if i > 105 or j > 105: 
                        continue

                    mcss1, mcss2, valid_mcss_pairs = find_mcss_matches(p1, self.sm1, p2, self.sm2)
                    if len(valid_mcss_pairs) == 0:
                        f.write('ERROR\n')
                        print 'no mcss found for poses', self.l1, self.l2, self.st
                        return

                    rmsd = 10000
                    for m1, m2 in valid_mcss_pairs:
                        nrmsd = calculate_in_place_rmsd(p1, mcss1[m1], p2, mcss2[m2], use_symmetry=True)
                        rmsd = min(rmsd, nrmsd)

                    f.write('{},{},{}\n'.format(i, j, rmsd))

if __name__ == '__main__':
    lig1 = sys.argv[1]
    lig2 = sys.argv[2]
    struct = sys.argv[3]
    glide_dir = 'docking/' + sys.argv[4]
    mcss_dir = '../../ligands/mcss/' + sys.argv[5]
    MCSS(lig1, lig2, struct, glide_dir, mcss_dir)






import os
import sys

def load_mcss(lig1, lig2):
    tr = {lig1:[],lig2:[]}
    m = None
    try:
        with open('{}-{}.csv'.format(lig1, lig2)) as f:
            for line in f:
                line = line.strip().split(',')
                if line[0] == 'SMILES': continue
                assert line[1] in tr
                lig,msize,smarts = line[1], int(line[5]), line[-1]
                tr[lig].append(smarts)
                if m is None: m = msize
                assert msize == m
    except Exception as e:
        print e
        print 'mcss file error'
    return m,tr[lig1],tr[lig2]

def delete_extraneous_bonds(st1, st2):
    
    if len(st1.bond) + 1 == len(st2.bond): # st2 has the extra bond
        st1, st2 = st2, st1

    all_b = {}
    for i,b in enumerate(st1.bond):
        all_b[(i,b.order)] = b 

    for (i,o),b in all_b.items():
        a1,a2 = b.atom1, b.atom2
        st1.deleteBond(a1,a2)
        if st1.isEquivalent(st2,False):
            return (st1,st2)
        st1.addBond(a1,a2,o)
    
    return None

def match(st1, st2):
    if abs(len(st1.bond) - len(st2.bond)) == 1:
        return delete_extraneous_bonds(st1, st2)
    if st1.isEquivalent(st2): return (st1,st2)
    return None

def find_mcss_matches(st1, smarts1, st2, smarts2):
    valid_pairs = []
    for sm1 in smarts1:
        for sm2 in smarts2:

            mcss1 = evaluate_smarts(st1, sm1, unique_sets=True)
            mcss2 = evaluate_smarts(st2, sm2, unique_sets=True)

            for m1 in mcss1:
                for m2 in mcss2:
                    equiv = match(st1.extract(m1), st2.extract(m2))
                    if equiv is not None: 
                        valid_pairs.append(equiv)

    return valid_pairs

ref_path = '../../../ligands/prepared_ligands/{}/{}.mae'
pv_path = '../../../docking/{}/{}-to-{}/{}-to-{}_pv.maegz'
outf = '{}-{}-to-{}-{}.csv'

class MCSS:
    def __init__(self, l1, l2, st=None, gdir=None, type_file=None):
        self.l1 = l1
        self.l2 = l2
        self.name = '{}-{}'.format(l1,l2)
        self.st = st
        self.gdir = gdir
        self.type_file = type_file

    def load_ref(self):
        ref1 = StructureReader(ref_path.format(self.l1, self.l1)).next()
        ref2 = StructureReader(ref_path.format(self.l2, self.l2)).next()
        return ref1, ref2

    def proc_ref(self):
        s3, sm1, sm2 = load_mcss(self.l1, self.l2)
        ref1, ref2 = self.load_ref()
        valid_mcss_pairs = find_mcss_matches(ref1, sm1, ref2, sm2)

        stwr = StructureWriter('{}-{}_mcss.mae'.format(self.l1, self.l2))
        for m1,m2 in valid_mcss_pairs:
            stwr.append(m1)
            stwr.append(m2)
        if len(valid_mcss_pairs) == 0:
            with open('log','a') as f:
                f.write('ERROR: reference validation')
            for ref, sm_list in {ref1:sm1,ref2:sm2}.items():
                for sm in sm_list: 
                    try: stwr.append(ref.extract(evaluate_smarts(ref, sm, unique_sets=True)[0]))
                    except:
                        with open('log','a') as f:
                            f.write('ERROR: evaluating smarts')
        stwr.close()

    def proc_pv(self):
        s3, sm1, sm2 = load_mcss(self.l1, self.l2)

        pv1 = list(StructureReader(pv_path.format(self.gdir, self.l1, self.st, self.l1, self.st)))[1:]
        pv2 = list(StructureReader(pv_path.format(self.gdir, self.l2, self.st, self.l2, self.st)))[1:]

        s1 = len([a for a in pv1[0].atom if a.element != 'H'])
        s2 = len([a for a in pv2[0].atom if a.element != 'H'])

        with open(outf.format(self.l1,self.l2,self.st,self.gdir), 'w') as f:
            f.write('{},{},{}\n'.format(s1, s2, s3))
            if s3 == 0: return

            for i, p1 in enumerate(pv1):
                for j, p2 in enumerate(pv2):
                    if i > 105 or j > 105: 
                        continue

                    valid_mcss_pairs = find_mcss_matches(p1, sm1, p2, sm2)
                    if len(valid_mcss_pairs) == 0:
                        f.write('ERROR\n')
                        print 'no mcss found for poses', self.l1, self.l2, self.st
                        return

                    rmsd = 10000
                    for m1, m2 in valid_mcss_pairs:
                        nrmsd = calculate_in_place_rmsd(m1, [a.index for a in m1.atom], 
                                                        m2, [a.index for a in m2.atom],  use_symmetry=True)
                        rmsd = min(rmsd, nrmsd)

                    f.write('{},{},{}\n'.format(i, j, rmsd))

if __name__ == '__main__':
    from schrodinger.structure import StructureReader, StructureWriter
    from schrodinger.structutils.rmsd import calculate_in_place_rmsd
    from schrodinger.structutils.analyze import evaluate_smarts

    mode = sys.argv[1]
    lig1 = sys.argv[2]
    lig2 = sys.argv[3]
    
    if mode == 'RMSD':
        mcss = MCSS(lig1, lig2, st=sys.argv[4], gdir=sys.argv[5])
        mcss.proc_pv()
    elif mode == 'REF': # for debugging
        mcss = MCSS(lig1, lig2)
        mcss.proc_ref()



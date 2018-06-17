import os
import sys
 
def load_mcss_csv(mpath):
    lig1 = mpath.split('-')[0]
    lig2 = mpath.split('-')[1]
    tr = {lig1:[],lig2:[]}
    m = None
    try:
        with open(mpath) as f:
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
        print 'mcss csv file error'
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
    if st1.isEquivalent(st2,False): return (st1,st2)
    return None

def unproc_st(st):
    st2 = st.copy()
    for b in st2.bond: b.order = 1
    for a in st2.atom: a.formal_charge = 0
    st2.retype()
    return st2

def find_mcss_matches(st1, smarts1, st2, smarts2, unproc=True):
    all_pairs = []
    for sm1 in smarts1:
        for sm2 in smarts2:

            mcss1 = evaluate_smarts(st1, sm1, unique_sets=True)
            mcss2 = evaluate_smarts(st2, sm2, unique_sets=True)

            if unproc:
                list1 = [unproc_st(st1.extract(m)) for m in mcss1]
                list2 = [unproc_st(st2.extract(m)) for m in mcss2]
            else:
                list1 = [st1.extract(m) for m in mcss1]
                list2 = [st2.extract(m) for m in mcss2]

            for m1 in list1:
                for m2 in list2:
                    attempt_match = match(m1,m2)
                    if attempt_match is not None:
                        all_pairs.append((m1,m2))

    return all_pairs

#ref_path = '../../../structures/ligands/{}.mae'
ref_path = '../../../ligands/prepared_ligands/{}/{}.mae'
pv_path = '../../../docking/{}/{}-to-{}/{}-to-{}_pv.maegz'
outf = '{}-{}-{}.csv'

class WritePairMCSS:
    def __init__(self, mpath, st=None, gdir=None):
        self.l1 = mpath.split('-')[0]
        self.l2 = mpath.split('-')[1]
        self.mpath = mpath
        self.st = st
        self.gdir = gdir

    def load_ref(self):
        ref1 = StructureReader(ref_path.format(self.l1, self.l1)).next()
        ref2 = StructureReader(ref_path.format(self.l2, self.l2)).next()
        return ref1, ref2

    def write_size_file(self):
        s3,sm1,sm2 = load_mcss_csv(self.mpath)
        ref1,ref2 = self.load_ref()
        valid_mcss_pairs = find_mcss_matches(ref1, sm1, ref2, sm2)

        s1 = len([a for a in ref1.atom if a.element != 'H'])
        s2 = len([a for a in ref2.atom if a.element != 'H'])

        with open('{}.size'.format(self.mpath.split('.')[0]),'w') as f:
            f.write('{} matches\n'.format(len(valid_mcss_pairs)))
            for smarts in sm1:
                f.write('{},{},{},{}\n'.format(self.l1,s1,s3,smarts))
            for smarts in sm2:
                f.write('{},{},{},{}\n'.format(self.l2,s2,s3,smarts))

    def proc_ref(self):
        s3, sm1, sm2 = load_mcss_csv(self.mpath)
        ref1, ref2 = self.load_ref()
        valid_mcss_pairs = find_mcss_matches(ref1, sm1, ref2, sm2)

        stwr = StructureWriter('{}.mae'.format(self.mpath.split('.')[0]))
        stwr.append(ref1)
        stwr.append(ref2)
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

    def proc_debug(self):
        s3, sm1, sm2 = load_mcss_csv(self.l1, self.l2, self.tf)
        print s3, sm1, sm2
        pv1 = list(StructureReader(pv_path.format(self.gdir, self.l1, self.st, self.l1, self.st)))[1:]
        pv2 = list(StructureReader(pv_path.format(self.gdir, self.l2, self.st, self.l2, self.st)))[1:]
        print 'hi'
        stwr = StructureWriter('debug.mae')
        for i, p1 in enumerate(pv1):
            for j, p2 in enumerate(pv2):
                if i > 3 or j > 3: 
                    continue
                print i,j
                all_mcss_pairs = find_mcss_matches(p1, sm1, p2, sm2)
                print len(all_mcss_pairs)
                stwr.append(p1)
                stwr.append(p2)
                rmsd = 10000
                count = 0
                for m1, m2 in all_mcss_pairs:
                    count += 1
                    stwr.append(m1)
                    stwr.append(m2)
                    try: 
                        conf_rmsd = ConformerRmsd(m1, m2)#.calculate()
                        conf_rmsd.use_heavy_atom_graph = True
                        conf_rmsd = conf_rmsd.calculate()
                    except: pass
                    print count, conf_rmsd
                    rmsd = min(rmsd, conf_rmsd)
        stwr.close()


    def write_rmsd_file(self):
        s3, sm1, sm2 = load_mcss_csv(self.mpath)

        pv1 = list(StructureReader(pv_path.format(self.gdir, self.l1, self.st, self.l1, self.st)))[1:]
        pv2 = list(StructureReader(pv_path.format(self.gdir, self.l2, self.st, self.l2, self.st)))[1:]

        s1 = len([a for a in pv1[0].atom if a.element != 'H'])
        s2 = len([a for a in pv2[0].atom if a.element != 'H'])

        with open(outf.format(self.mpath.split('.')[0],self.st,self.gdir), 'w') as f:
            f.write('{},{},{}\n'.format(s1, s2, s3))
            if s3 <= 8: return # or s3*2 < min(s1,s2): return

            for i, p1 in enumerate(pv1):
                for j, p2 in enumerate(pv2):
                    if i > 105 or j > 105: 
                        continue
                    
                    rmsd = 10000
                    for m1, m2 in find_mcss_matches(p1,sm1,p2,sm2):
                        try: 
                            conf_rmsd = ConformerRmsd(m1, m2)
                            conf_rmsd.use_heavy_atom_graph = True
                            rmsd = min(rmsd, conf_rmsd.calculate())
                        except: pass

                    if rmsd == 10000:
                        f.write('ERROR\n')
                        print i,j,'no mcss found',self.l1,self.l2,self.st
                        return

                    f.write('{},{},{}\n'.format(i, j, rmsd))

if __name__ == '__main__':
    from schrodinger.structure import StructureReader, StructureWriter
    from schrodinger.structutils.rmsd import ConformerRmsd
    from schrodinger.structutils.analyze import evaluate_smarts

    mode = sys.argv[1]
    mpath = sys.argv[2]
    if mode == 'RMSD':
        mcss = WritePairMCSS(mpath, st=sys.argv[3], gdir=sys.argv[4])
        mcss.write_rmsd_file()
    elif mode == 'SIZE':
        mcss = WritePairMCSS(mpath)
        mcss.write_size_file()
    elif mode == 'REF':
        mcss = WritePairMCSS(mpath)
        mcss.proc_ref()
    elif mode == 'DEBUG':
        mcss = WritePairMCSS(mpath, st='1E3G', gdir='glide12')
        mcss.proc_debug()


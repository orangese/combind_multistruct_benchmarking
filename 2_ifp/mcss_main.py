import os
import sys

class csvFile:
    def __init__(self, mpath):
        self.mpath = mpath
        self.l1 = mpath.split('-')[0]
        self.l2 = mpath.split('-')[1]
        self.msz = -1
        self.sm = {}
        self.alist = {}
        self.num_b = -1        
        self.load()
        
    def load(self):
        try:
            with open(self.mpath) as f:
                for line in f:
                    csv = line.strip().split(',')
                    qsv = line.strip().split('"')
                    if csv[0] == 'SMILES': continue
                    lig = csv[1]
                    assert lig in self.mpath
                    if lig not in self.sm: self.sm[lig] = []
                    if lig not in self.alist: self.alist[lig] = []
                    self.sm[lig].append(csv[-1])
                    if len(qsv[1]) > 0: self.alist[lig].append(tuple(sorted([int(a) for a in qsv[1].split(',')])))
                    if self.msz == -1: self.msz = int(csv[5])
                    if self.num_b == -1: self.num_b = int(csv[6])
                    assert int(csv[6]) == self.num_b
                    assert int(csv[5]) == self.msz#msize == m
        except Exception as e:
            print e
            print self.mpath
            print 'mcss csv file error'

def csv_equality(f1, f2, smarts=False, alist=False, size=False):
    match = True
    if size:
        match = match and f1.num_b == f2.num_b and f1.msz == f2.msz

    if smarts:
        match = match and sorted(f1.sm.keys()) == sorted(f2.sm.keys())
        for lig,sm1 in f1.sm.items():
            _sm1 = f2.sm[lig]
            match = match and sorted(sm1) == sorted(_sm1)
    
    if alist:
        match = match and sorted(f1.alist.keys()) == sorted(f2.alist.keys())
        for lig,alist1 in f1.alist.items():
            _alist1 = f2.alist[lig]
            match = match and sorted(alist1) == sorted(_alist1)

    return match

def delete_extraneous_bonds(st1, st2):    
    if len(st1.bond) + 1 == len(st2.bond): # st2 has the extra bond
        st1, st2 = st2, st1

    all_b = {(i,b.order):b for i,b in enumerate(st1.bond)}
    for (i,o),b in all_b.items():
        a1,a2 = b.atom1, b.atom2
        st1.deleteBond(a1,a2)
        if st1.isEquivalent(st2,False):
            return (st1,st2)
        st1.addBond(a1,a2,o)

def delete_bond_pair(st1, st2):
    all_b1 = {(i,b.order):b for i,b in enumerate(st1.bond)}
    all_b2 = {(i,b.order):b for i,b in enumerate(st2.bond)}
    for (i,o1),b1 in all_b1.items():
        a1,a2 = b1.atom1,b1.atom2
        st1.deleteBond(a1,a2)
        for (j,o2),b2 in all_b2.items():
            a3,a4 = b2.atom1,b2.atom2
            st2.deleteBond(a3,a4)
            if st1.isEquivalent(st2,False):
                return (st1,st2)
            st2.addBond(a3,a4,o2)
        st1.addBond(a1,a2,o1)

def match(st1, st2, num_b):
    if st1.isEquivalent(st2,False): return (st1,st2)
    if abs(len(st1.bond) - len(st2.bond)) == 1:
        return delete_extraneous_bonds(st1, st2)
    if len(st1.bond) == len(st2.bond) and len(st1.bond) == num_b + 1:
        return delete_bond_pair(st1, st2)

def unproc_st(st):
    st2 = st.copy()
    for b in st2.bond: b.order = 1
    for a in st2.atom: a.formal_charge = 0
    st2.retype()
    return st2

def find_mcss_matches(st1, st2, csvf, unproc=True):
    all_pairs = []
    
    for sm1 in csvf.sm[csvf.l1]:
        for sm2 in csvf.sm[csvf.l2]:
            if '[11C]' in sm1 or '[11C]' in sm2: return []
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
                    attempt_match = match(m1,m2,csvf.num_b)
                    if attempt_match is not None:
                        all_pairs.append((m1,m2))
    return all_pairs

ref_path2 = '../../../structures/ligands/{}.mae'
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
        self.csv_f = csvFile(mpath)

    def load_ref(self):#,alt=False):
        #if alt: ref_path = ref_path2
        ref1 = StructureReader(ref_path.format(self.l1, self.l1)).next()
        ref2 = StructureReader(ref_path.format(self.l2, self.l2)).next()
        return ref1, ref2

    def write_size_file(self):
        other_mcss_files = [f.split('.')[0]+'.csv' for f in os.listdir('.') if f[-4:] == 'size']
        for o in sorted(other_mcss_files):
            if o == self.mpath: continue
            of = csvFile(o)
            if csv_equality(self.csv_f,of,size=True,alist=True):
                more_f = [f for f in os.listdir('.') if o.split('.')[0] in f.split('.')[0]]
                if True in [f.split('.')[-1] == 'size' for f in more_f]: 
                    os.system('cp {}.size {}.size'.format(o.split('.')[0],self.mpath.split('.')[0]))
                    for f in more_f:
                        if 'glide' in f:
                            other_rmsd_info = '-'.join(f.split('-')[-2:])
                            new_rmsd_f = '{}-{}'.format(self.mpath.split('.')[0],other_rmsd_info)
                            os.system('cp {} {}'.format(f,new_rmsd_f))
                    return

        ref1,ref2 = self.load_ref()
       
        s1 = len([a for a in ref1.atom if a.element != 'H'])
        s2 = len([a for a in ref2.atom if a.element != 'H'])
 
        if self.csv_f.msz*2 < min(s1,s2):
            num_m = -1
        else:
            num_m = len(find_mcss_matches(ref1, ref2, self.csv_f))

        with open('{}.size'.format(self.mpath.split('.')[0]),'w') as f:
            f.write('{} matches\n'.format(num_m))
            for smarts in self.csv_f.sm[self.l1]:
                f.write('{},{},{},{}\n'.format(self.l1,s1,self.csv_f.msz,smarts))
            for smarts in self.csv_f.sm[self.l2]:
                f.write('{},{},{},{}\n'.format(self.l2,s2,self.csv_f.msz,smarts))

    def proc_ref(self):
        ref1, ref2 = self.load_ref()#alt=True)
        valid_mcss_pairs = find_mcss_matches(ref1, ref2, self.csv_f)

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
        pv1 = list(StructureReader(pv_path.format(self.gdir, self.l1, self.st, self.l1, self.st)))[1:]
        pv2 = list(StructureReader(pv_path.format(self.gdir, self.l2, self.st, self.l2, self.st)))[1:]
        print 'hi'
        stwr = StructureWriter('debug.mae')
        for i, p1 in enumerate(pv1):
            for j, p2 in enumerate(pv2):
                if i > 3 or j > 3: 
                    continue
                print i,j
                all_mcss_pairs = find_mcss_matches(p1, sm1, p2, sm2, num_b)
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
        
        outpth = outf.format(self.mpath.split('.')[0],self.st,self.gdir)

        pv1 = list(StructureReader(pv_path.format(self.gdir, self.l1, self.st, self.l1, self.st)))[1:]
        pv2 = list(StructureReader(pv_path.format(self.gdir, self.l2, self.st, self.l2, self.st)))[1:]

        s1 = len([a for a in pv1[0].atom if a.element != 'H'])
        s2 = len([a for a in pv2[0].atom if a.element != 'H'])

        with open(outpth, 'w') as f:
            f.write('{},{},{}\n'.format(s1, s2, self.csv_f.msz))
            if self.csv_f.msz*2 < min(s1,s2): return

            for i, p1 in enumerate(pv1):
                for j, p2 in enumerate(pv2):
                    if i > 105 or j > 105: 
                        continue
                    
                    rmsd = 10000
                    for m1, m2 in find_mcss_matches(p1,p2,self.csv_f):
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
    from schrodinger.structutils.analyze import evaluate_smarts
    from schrodinger.structutils.rmsd import ConformerRmsd

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


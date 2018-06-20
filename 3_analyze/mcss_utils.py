import os
import sys

class MCSS:
    def __init__(self, sp, st, root, num_poses=None):
        self.sp = sp
        self.st = st
        self.root = root
        self.num_poses = num_poses
        self.all_pairs = {}

    def get_pair(self, l1, l2):
        if (l1,l2) in self.all_pairs: return self.all_pairs[(l1,l2)]
        if (l2,l1) in self.all_pairs: return self.all_pairs[(l2,l1)]

    def load_mcss(self, lset1, lset2, rmsd=False):
        for l1 in lset1:
            for l2 in lset2:
                if l1 == l2: continue
                if (l1,l2) not in self.all_pairs:
                    if os.path.exists('{}/mcss/{}/{}-{}'.format(self.root, self.sp['mcss'], l1, l2)):
                        self.all_pairs[(l1,l2)] = PairMCSS(l1,l2,self.st,self.sp,self.root)
                        self.all_pairs[(l1,l2)].load(rmsd,self.num_poses)
                    elif os.path.exists('{}/mcss/{}/{}-{}'.format(self.root, self.sp['mcss'], l2, l1)):
                        self.all_pairs[(l2,l1)] = PairMCSS(l2,l1,self.st,self.sp,self.root)
                        self.all_pairs[(l2,l1)].load(rmsd,self.num_poses)

    def sort_by_mcss(self, q, lset, sort_prop=lambda x: x.m_sz):
        tr = {}
        for p in [p for lp,p in self.all_pairs.items() if q in lp]:# and p.num_matches > 0]:
            if p.l1 == q and p.l2 in lset: tr[p.l2] = p
            if p.l2 == q and p.l1 in lset: tr[p.l1] = p
        return sorted(tr.keys(),key=lambda x: -sort_prop(tr[x]))

    def get_rmsd(self, l1, p1, l2, p2):
        if (l1,l2) in self.all_pairs: return self.all_pairs[(l1,l2)].rmsds.get((p1,p2),None)
        if (l2,l1) in self.all_pairs: return self.all_pairs[(l2,l1)].rmsds.get((p2,p1),None)
        return None

class PairMCSS:
    def __init__(self, l1, l2, st, sp, root):
        self.l1 = l1
        self.l2 = l2
        self.st = st
        self.sp = sp
        self.root = root
        self.name = '{}-{}'.format(l1,l2)

        self.m_sz = 0
        self.l_sz = {l1:0,l2:0}
        self.num_matches = 0
        self.smarts = {l1:[],l2:[]}    

        self.rmsds = {}
        
    def get_mcss_path(self,add_dir=False,add_st=False,ext='csv'):
        pth = '{}-{}'.format(self.name,self.sp['mcss_type'])
        if add_st: pth = '{}-{}-{}'.format(pth,self.st,self.sp['docking'])
        if add_dir: pth = '{}/mcss/{}/{}/{}'.format(self.root,self.sp['mcss'],self.name,pth)
        return '{}.{}'.format(pth,ext)

    def load(self,rmsd,num_poses):
        self.load_mcss_size()
        if rmsd and self.m_sz >= 16 and self.num_matches > 0: 
            self.load_rmsds(num_poses)

    def load_mcss_size(self):
        pth = self.get_mcss_path(add_dir=True,ext='size')
        if not os.path.exists(pth):
            #print 'no size file',pth
            #os.system('rm -f {}'.format(self.get_mcss_path(add_dir=True,add_st=True)))
            return False
        try:
            with open(pth) as f:
                for i,line in enumerate(f):
                    if i == 0:
                        self.num_matches = int(line.split(' ')[0])
                        #if self.num_matches == 0:
                            #os.system('rm -f {}'.format(pth))
                            #print 'warning: no matches found for ', self.l1,self.l2,self.sp['mcss_type']
                            #break
                        continue
                    lig,lsize,msize,smarts = line.strip().split(',')
                    self.m_sz = int(msize)
                    self.l_sz[lig] = int(lsize)
                    self.smarts[lig].append(smarts)
        except Exception as e:
            print e
            print 'mcss size file error'
            return False
        return True

    def load_rmsds(self, num_poses):
        if len(self.rmsds) > 0: return
        pth = self.get_mcss_path(add_dir=True,add_st=True)
        if os.path.exists(pth) and os.stat(pth).st_size == 0:
            print 'deleting empty file',pth
            return
        if not os.path.exists(pth): return#continue
        try:
            with open(pth) as f:
                for line_count, line in enumerate(f):
                    line = line.strip().split(',')
                    if line_count == 0:
                        s1, s2, s3 = [float(i) for i in line]
                        assert s1 == self.l_sz[self.l1]
                        assert s2 == self.l_sz[self.l2]
                        assert s3 == self.m_sz
                    if 'ERROR' in line:
                        print 'ERROR',self.l1,self.l2,s1, s2, s3
                        break
                    p1,p2,rmsd = int(line[0]), int(line[1]), float(line[2])
                    self.rmsds[(p1,p2)] = rmsd
                else:
                    if line_count < 2:
                        print line_count,'too small'
                    elif p1+1 < min(num_poses[self.l1], 100) or p2+1 < min(num_poses[self.l2], 100) or rmsd == float(10000):
                        print line_count, 'hmmmm', pth, p1, p2, num_poses[self.l1], num_poses[self.l2], rmsd
                        #os.system('rm {}'.format(fpath))
        except Exception as e:
            print 'mcss parse error', pth
            print e
            print line
            #os.system('rm {}'.format(fpath))


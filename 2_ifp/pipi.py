from schrodinger.structutils.measure import measure_distance
from schrodinger.structutils.analyze import center_of_mass

class PiPi:
    def __init__(self, r1, i1, r2, i2): # type = structure
        self.r1 = r1#[r for r in r1.ring][0] # type = ring
        self.r2 = r2#[r for r in r2.ring][0]

        self.i1 = i1 # protein ring unique id
        self.i2 = i2 # ligand ring unique id

        self.dist = self.get_dist()
        self.lig_ring = self.r2

    def get_dist(self):
        st1 = self.r1.extractStructure(copy_props=True)
        st2 = self.r2.extractStructure(copy_props=True)
        return measure_distance(center_of_mass(st1), center_of_mass(st2))

    def score(self):
        if self.dist <= 6: return 1
        elif self.dist <= 8: return (8.0-self.dist)/2.0
        else: return 0

    def newscore(self):
        if self.dist <= 4.5: return 1
        elif self.dist <= 6: return (6.0-self.dist)/1.5
        return 0

    def __str__(self):
        str1= '\nPiPi:\nscore: ' + str(self.score())+','+ str(self.newscore())+'\n'
        str2= '+ring 1: \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.r1.atom), self.r1.isAromatic(), self.r1.isHeteroaromatic())
        str3= '+ring 2: \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.r2.atom), self.r2.isAromatic(), self.r2.isHeteroaromatic())
        return str1+str2+str3

class PiPi_Container:
    def __init__(self, lig_st, sub_st_map, ind):
        self.lig_st = lig_st
        self.sub_st_map = sub_st_map
        self.ind = ind
        #self.lig_aro = [ri.extractStructure(copy_props=True) for ri in lig_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        self.lig_aro = [ri for ri in lig_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        #for i, ri in enumerate(lig_st.ring):
        #    print i, ri.isAromatic(), ri.isHeteroaromatic()
        #    print sorted([(a.index, a.element) for a in ri.atom]), self.assign_ss_to_ring(ri)
        #    b = ri.extractStructure(copy_props=True)
        #    print sorted([(a.index, a.element) for a in b.atom]), self.assign_ss_to_ring(b)
        self.all_pipi = {}

    def add_residue(self, resnum, res_st):
        self.all_pipi[resnum] = []
        res_aro = [ri for ri in res_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        for i1, r1 in enumerate(res_aro):
            for i2, r2 in enumerate(self.lig_aro):
                self.all_pipi[resnum].append(PiPi(r1, (i1,resnum), r2, i2))

    def filter_int(self):
        filtered_pipi = {}
        for r in self.all_pipi:
            for pipi in self.all_pipi[r]:
                filtered_pipi[(r, pipi.i1, pipi.i2)] = pipi
        ranked_pipi = sorted(filtered_pipi.keys(), key=lambda x: -filtered_pipi[x].score() - filtered_pipi[x].newscore())
        self.all_pipi = {}
        used_i1 = set()
        used_i2 = set()
        for pipi_key in ranked_pipi:
            r, i1, i2 = pipi_key
            if i1 not in used_i1 and i2 not in used_i2 and filtered_pipi[pipi_key].score() > 0:
                if r not in self.all_pipi:
                    self.all_pipi[r] = []
                self.all_pipi[r].append(filtered_pipi[pipi_key])
                used_i1.add(i1)
                used_i2.add(i2) 

    def assign_ss_to_ring(self, ring):
        all_ss = []
        for a in ring.atom:
            if a.element == 'H': continue
            all_ss.extend(self.sub_st_map.get(a.index, ['']))
        #print 'ring', self.lig_st._getTitle(), all_ss
        return str(sorted(all_ss, key=lambda x:-len([i for i in all_ss if i == x]))[0])

    def score(self):
        all_scores = {}
        for r, pipi_list in self.all_pipi.items():
            for p in pipi_list:
                majority_ss = ''# '{},{}'.format(p.i1, p.i2)#majority_ss = self.assign_ss_to_ring(p.lig_ring)
                key1 = (self.ind[0], r, majority_ss)
                key2 = (self.ind[1], r, majority_ss)
                all_scores[key1] = all_scores.get(key1, 0) + p.score()
                all_scores[key2] = all_scores.get(key2, 0) + p.newscore()
        return all_scores

        #return {
        #    r : [
        #        sum([pipi.score() for pipi in self.all_pipi[r]]),
        #        sum([pipi.newscore() for pipi in self.all_pipi[r]])
        #    ] for r in self.all_pipi
        #}

    def __str__(self):
        return '\nAll PiPi:\n' + '\n'.join([str(r) + '\n'.join([str(i) for i in self.all_pipi[r]]) for r in sorted(self.all_pipi.keys())])


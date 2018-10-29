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
        if self.dist <= params['pipi_dist_opt']:
            return 1
        elif self.dist <= params['pipi_dist_cut']:
            return ((params['pipi_dist_cut'] - self.dist)
                    / (params['pipi_dist_cut'] - params['pipi_dist_opt']))
        else:
            return 0

    def __str__(self):
        str1= '\nPiPi:\nscore: ' + str(self.score())+'\n'
        str2= '+ring 1: \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.r1.atom), self.r1.isAromatic(), self.r1.isHeteroaromatic())
        str3= '+ring 2: \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.r2.atom), self.r2.isAromatic(), self.r2.isHeteroaromatic())
        return str1+str2+str3

class PiPi_Container:
    def __init__(self, lig, ind):
        self.lig = lig
        self.ind = ind
        self.all_pipi = {}

    def add_residue(self, resnum, res):
        self.all_pipi[resnum] = []
        for i1, r1 in enumerate(res.aro):
            for i2, r2 in enumerate(self.lig.aro):
                self.all_pipi[resnum].append(PiPi(r1, (i1,resnum), r2, i2))

    def filter_int(self):
        filtered_pipi = {}
        for r in self.all_pipi:
            for pipi in self.all_pipi[r]:
                filtered_pipi[(r, pipi.i1, pipi.i2)] = pipi
        ranked_pipi = sorted(filtered_pipi.keys(), key=lambda x: -filtered_pipi[x].score())
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

    def score(self):
        all_scores = {}
        for r, pipi_list in self.all_pipi.items():
            for p in pipi_list:
                key1 = (self.ind[0], r, '')
                all_scores[key1] = all_scores.get(key1, 0) + p.score()
        return all_scores

    def __str__(self):
        return '\nAll PiPi:\n' + '\n'.join([str(r) + '\n'.join([str(i) for i in self.all_pipi[r]]) for r in sorted(self.all_pipi.keys())])


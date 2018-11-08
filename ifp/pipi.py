from schrodinger.structutils.measure import measure_distance
from schrodinger.structutils.analyze import center_of_mass
from shared_paths import shared_paths
params = shared_paths['ifp']

class PiPi_Container:
    def __init__(self, lig, ind):
        self.lig = lig
        self.ind = ind
        self.all_pipi = {}

    def add_residue(self, resnum, res):
        self.all_pipi[resnum] = []
        for i1, r1 in enumerate(res.aro):
            for i2, r2 in enumerate(self.lig.aro):
                self.all_pipi[resnum].append(PiPi(r1, (i1, resnum), r2, i2))

    def filter_int(self):
        pass

    def score(self):
        all_scores = {}
        for r, pipi_list in self.all_pipi.items():
            for p in pipi_list:
                key = (self.ind[0], r, '')
                if key not in all_scores: all_scores[key] = 0
                all_scores[key] += p.score()
        return all_scores

    def raw(self):
        all_raw = {}
        for r, pipi_list in self.all_pipi.items():
            for p in pipi_list:
                if not p.score(): continue
                key = (self.ind[0], r, '')
                if key not in all_raw: all_raw[key] = []
                all_raw[key] += [p.dist]
        return all_raw

class PiPi:
    """
    Assesses the potential for Pi-Pi interaction between the substructures
    r1 and r2 identified by i1 and i2 respectively. It is assumed that
    r1 and r2 each contain 1 or more aromatic rings. The minimum distance
    between the center of mass of any aromatic ring in r1 and any aromatic
    ring in r2 is used to determine the distance between the two.

    The input rings are expected to be generated through the AtomGroup.fuse_rings
    method.
    """
    def __init__(self, r1, i1, r2, i2): # type = structure
        self.r1 = r1
        self.r2 = r2

        self.i1 = i1 # protein ring unique id
        self.i2 = i2 # ligand ring unique id

        self.dist = self.get_dist()

    def get_dist(self):
        dist = float('inf')
        for r1 in self.r1.ring:
            for r2 in self.r2.ring:
                st1 = r1.extractStructure(copy_props=True)
                st2 = r2.extractStructure(copy_props=True)
                dist = min(dist, measure_distance(center_of_mass(st1), center_of_mass(st2)))
        assert dist != float('inf')
        return dist

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

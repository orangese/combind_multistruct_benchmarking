from schrodinger.structutils.measure import measure_distance
from schrodinger.structutils.analyze import center_of_mass

class PiPi:
    def __init__(self, r1, r2): # type = structure
        self.r1 = [r for r in r1.ring][0] # type = ring
        self.r2 = [r for r in r2.ring][0]
        self.dist = measure_distance(center_of_mass(r1), center_of_mass(r2))#self.c1.dist_to(self.c2)

    def score(self):
        if self.dist <= 6: return 1
        elif self.dist <= 8: return (8 - self.dist)/2.0
        else: return 0

    def __str__(self):
        str1= '\nPiPi:\nscore: ' + str(self.score())+'\n'
        str2= '+ring 1: \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.r1.atom), self.r1.isAromatic(), self.r1.isHeteroaromatic())
        str3= '+ring 2: \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.r2.atom), self.r2.isAromatic(), self.r2.isHeteroaromatic())
        return str1+str2+str3

class PiPi_Container:
    def __init__(self, lig_st):
        self.lig_st = lig_st
        self.lig_aro = [ri.extractStructure() for ri in lig_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]

        self.all_pipi = {}

    def add_residue(self, resnum, res_st):
        self.all_pipi[resnum] = []
        res_aro = [ri.extractStructure() for ri in res_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        for r1 in res_aro:
            for r2 in self.lig_aro:
                self.all_pipi[resnum].append(PiPi(r1, r2))

    def filter_int(self):
        self.all_pipi = {r:[pipi for pipi in self.all_pipi[r] if pipi.score() > 0.05] for r in self.all_pipi}
        self.all_pipi = {r:self.all_pipi[r] for r in self.all_pipi if len(self.all_pipi[r]) > 0}    

    def score(self):
        return {
            r : [
                sum([pipi.score() for pipi in self.all_pipi[r]])
            ] for r in self.all_pipi
        }

    def __str__(self):
        return '\nAll PiPi:\n' + '\n'.join([str(r) + '\n'.join([str(i) for i in self.all_pipi[r]]) for r in sorted(self.all_pipi.keys())])


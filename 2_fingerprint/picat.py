from schrodinger.structutils.measure import measure_distance
from schrodinger.structutils.analyze import center_of_mass, find_rings

class PiCat:
    def __init__(self, ri, cation, resIsPi):
        self.ring = [r for r in ri.ring][0] # type = ring
        self.cation = cation
        self.dist = measure_distance(center_of_mass(ri), cation.xyz)
        self.resIsPi = resIsPi

    def score(self):
        if self.dist <= 6: return 1
        elif self.dist <= 8: return (8 - self.dist)/2
        else: return 0

    def __str__(self):
        str1= 'PiCat:  \nscore: ' + str(self.score()) + '\n residue is Pi: {}\n'.format(self.resIsPi)
        str2= '+ring : \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.ring.atom), self.ring.isAromatic(), self.ring.isHeteroaromatic())
        str3= '+cation: \n++{} {}\n'.format(self.cation.element, self.cation.index)
        return str1+str2+str3


class PiCat_Container:
    def __init__(self, lig_st):
        self.lig_st = lig_st
        self.lig_aro = [ri.extractStructure() for ri in lig_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        self.lig_cat = [a for a in lig_st.atom if a.formal_charge > 0]
        self.all_picat = {}

    def add_residue(self, resnum, res_st):
        self.all_picat[resnum] = []
        res_aro = [ri.extractStructure() for ri in res_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        res_cat = [a for a in res_st.atom if a.formal_charge > 0]
        
        for aro in res_aro:
            for cat in self.lig_cat:
                self.all_picat[resnum].append(PiCat(aro, cat, True))

        for cat in res_cat:
            for aro in self.lig_aro:
                self.all_picat[resnum].append(PiCat(aro, cat, False))

    def filter_int(self):
        self.all_picat = {r:[picat for picat in self.all_picat[r] if picat.score() > 0.05] for r in self.all_picat}
        self.all_picat = {r:self.all_picat[r] for r in self.all_picat if len(self.all_picat[r]) > 0}    

    def score(self):
        return {
            r : [
                sum([picat.score() for picat in self.all_picat[r] if picat.resIsPi]),
                sum([picat.score() for picat in self.all_picat[r] if not picat.resIsPi])
            ] for r in self.all_picat
        }

    def __str__(self):
        return 'Cation-Pi: '+''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_picat[r]]) for r in sorted(self.all_picat.keys())])


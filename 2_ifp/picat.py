from schrodinger.structutils.measure import measure_distance
from schrodinger.structutils.analyze import center_of_mass, find_rings

class PiCat:
    def __init__(self, ri, ri_ind, cation, cat_ind, resIsPi):
        self.ring = ri#[r for r in ri.ring][0] # type = ring
        self.cation = cation
        self.ring_ind = ri_ind
        self.cation_ind = cat_ind
        st1 = ri.extractStructure(copy_props=True)
        self.dist = measure_distance(center_of_mass(st1), cation.xyz)
        self.resIsPi = resIsPi

    def score(self):
        if self.dist <= 4.5: return 1
        elif self.dist <= 6: return (6 - self.dist)/1.5
        else: return 0

    def __str__(self):
        str1= 'PiCat:  \nscore: ' + str(self.score()) + '\n residue is Pi: {}\n'.format(self.resIsPi)
        str2= '+ring : \n++{} atoms\n++{} aromatic\n++{} heteroaromatic\n'.format(len(self.ring.atom), self.ring.isAromatic(), self.ring.isHeteroaromatic())
        str3= '+cation: \n++{} {}\n'.format(self.cation.element, self.cation.index)
        return str1+str2+str3


class PiCat_Container:
    def __init__(self, lig, ind):
        self.lig = lig
        self.ind = ind
        self.all_picat = {}

    def add_residue(self, resnum, res):
        self.all_picat[resnum] = []
        
        for i,aro in enumerate(res.aro):
            for j,cat in enumerate(self.lig.cat):
                self.all_picat[resnum].append(PiCat(aro, (i,resnum), cat, j, True))

        for j,cat in enumerate(res.cat):
            for i,aro in enumerate(self.lig.aro):
                self.all_picat[resnum].append(PiCat(aro, i, cat, (j,resnum), False))

    def filter_int(self):
        filtered_picat = {}
        for r in self.all_picat:
            for picat in self.all_picat[r]:
                filtered_picat[(r, picat.ring_ind, picat.cation_ind)] = picat
        ranked_picat = sorted(filtered_picat.keys(), key=lambda x: -filtered_picat[x].score())
        self.all_picat = {}
        used_i1 = set()
        used_i2 = set()
        for picat_key in ranked_picat:
            r, i1, i2 = picat_key
            if i1 not in used_i1 and i2 not in used_i2 and filtered_picat[picat_key].score() > 0:
                if r not in self.all_picat:
                    self.all_picat[r] = []
                self.all_picat[r].append(filtered_picat[picat_key])
                used_i1.add(i1)
                used_i2.add(i2)

    def score(self):
        all_scores = {}
        for r, pc_list in self.all_picat.items():#all_res():
            for pc in pc_list:
                if pc.resIsPi:
                    key = (self.ind[0], r, '')#lig_sub_st)
                else:
                    key = (self.ind[1], r, '')#lig_sub_st)
                all_scores[key] = all_scores.get(key, 0) + pc.score()
        return all_scores

    def __str__(self):
        return 'Cation-Pi: '+''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_picat[r]]) for r in sorted(self.all_picat.keys())])


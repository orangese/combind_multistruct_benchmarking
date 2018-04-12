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
    def __init__(self, lig_st, sub_st_map, ind):
        self.lig_st = lig_st
        self.sub_st_map = sub_st_map
        self.ind = ind
        self.lig_aro = [ri for ri in lig_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        self.lig_cat = [a for a in lig_st.atom if a.formal_charge > 0]
        self.all_picat = {}

    def add_residue(self, resnum, res_st):
        self.all_picat[resnum] = []
        res_aro = [ri for ri in res_st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        res_cat = [a for a in res_st.atom if a.formal_charge > 0]
        
        for i,aro in enumerate(res_aro):
            for j,cat in enumerate(self.lig_cat):
                self.all_picat[resnum].append(PiCat(aro, (resnum,i), cat, ('lig',j), True))

        for j,cat in enumerate(res_cat):
            for i,aro in enumerate(self.lig_aro):
                self.all_picat[resnum].append(PiCat(aro, ('lig',i), cat, (resnum,j), False))

    def filter_int(self):
        picat_list = []
        for r,r_picat in self.all_picat.items():
            picat_list.extend(r_picat)
        picat_list.sort(key=lambda x: -x.score()) 
        self.all_picat = {}
        used_ri = set()
        used_ca = set()
        for picat in picat_list:
            i1, i2 = picat.ring_ind, picat.cation_ind
            r = i1[0] if picat.resIsPi else i2[0]
            if i1 not in used_ri and i2 not in used_ca and picat.score() > 0:
                if r not in self.all_picat:
                    self.all_picat[r] = []
                self.all_picat[r].append(picat)
                used_ri.add(i1)
                used_ca.add(i2)

    def assign_ss_to_ring(self, ring):
        all_ss = []
        for a in ring.atom:
            if a.element == 'H': continue
            all_ss.extend(self.sub_st_map.get(a.index, ['']))
        #print 'ring', all_ss
        return str(sorted(all_ss, key=lambda x:-len([i for i in all_ss if i == x]))[0])

    def score(self):
        all_scores = {}
        for r, pc_list in self.all_picat.items():#all_res():
            for pc in pc_list:
                if pc.resIsPi:
                    lig_sub_st = ','.join([str(j) for j in self.sub_st_map.get(pc.cation.index, '')])
                    key = (self.ind[0], r, lig_sub_st)
                else:
                    lig_sub_st = self.assign_ss_to_ring(pc.ring)
                    key = (self.ind[1], r, lig_sub_st)
                all_scores[key] = all_scores.get(key, 0) + pc.score()
        return all_scores
        #return {
        #    r : [
        #        sum([picat.score() for picat in self.all_picat[r] if picat.resIsPi]),
        #        sum([picat.score() for picat in self.all_picat[r] if not picat.resIsPi])
        #    ] for r in self.all_picat
        #}

    def __str__(self):
        return 'Cation-Pi: '+''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_picat[r]]) for r in sorted(self.all_picat.keys())])


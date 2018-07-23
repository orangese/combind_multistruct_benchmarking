import math
from schrodinger.structutils.measure import measure_distance, measure_bond_angle

class HBond_Container:
    def __init__(self, lig, ind):
        self.lig = lig
        self.ind = ind
        self.all_hdon = {}
        self.all_hacc = {}

    def add_residue(self, resnum, res_st):
        self.all_hdon[resnum] = []
        self.all_hacc[resnum] = []

        for res_atom in res_st.hdon:
            for lig_atom in self.lig.hacc:
                self.all_hdon[resnum].extend([HBond(res_atom, lig_atom, n, resnum, True) for n in res_atom.bonded_atoms if n.element == 'H'])

        for res_atom in res_st.hacc:
            for lig_atom in self.lig.hdon:
                self.all_hacc[resnum].extend([HBond(lig_atom, res_atom, n, resnum, False) for n in lig_atom.bonded_atoms if n.element == 'H'])

    def filter_int(self):
        # enforces 1 hbond per h        
        unique_h = {}
        for r in set(list(self.all_hdon.keys()) + list(self.all_hacc.keys())):
            for hb in self.all_hdon.get(r, []) + self.all_hacc.get(r, []):
                if hb.h.index not in unique_h or unique_h[hb.h.index].score() < hb.score():
                    unique_h[hb.h.index] = hb

        self.all_hdon = {}
        self.all_hacc = {}
        for h, hb in unique_h.items():
            if hb.resIsHDonor:
                if hb.r_ind not in self.all_hdon: self.all_hdon[hb.r_ind] = []
                self.all_hdon[hb.r_ind].append(hb)
            if not hb.resIsHDonor:
                if hb.r_ind not in self.all_hacc: self.all_hacc[hb.r_ind] = []
                self.all_hacc[hb.r_ind].append(hb)

    def score(self):
        all_scores = {}
        for r, hb_list in self.all_hdon.items():#all_res():
            for hb in hb_list:
                key = (self.ind[0], r, '')
                all_scores[key] = all_scores.get(key, 0) + hb.score()
        for r, hb_list in self.all_hacc.items():
            for hb in hb_list:#self.all_hacc.get(r, []):
                key = (self.ind[1], r, '')
                all_scores[key] = all_scores.get(key, 0) + hb.score()
        return all_scores

class HBond:
    def __init__(self, donor, acceptor, h, r_ind, resIsHDonor): # D - H ... A - X
        self.d = donor # h donor, \in {N,O}
        self.a = acceptor # h acceptor, \in {N,O}
        self.h = h # covalently bound to the donor
        self.resIsHDonor = resIsHDonor

        self.r_ind = r_ind

        self.dist = measure_distance(h, acceptor) # angstroms
        self.DHA_angle = 180 - measure_bond_angle(self.d, self.h, self.a)#self.get_DHA_angle() # degrees

        self.charge_score = 1#1.2 if (self.d.formal_charge > 0 or self.a.formal_charge < 0) else 1

    def dist_score(self):
        if self.dist <= 2.5: return 1
        elif self.dist <= 3: return (3 - self.dist)/0.5
        else: return 0

    def angle_score(self):
        #if self.a.element == 'O':
        if self.DHA_angle <= 60: return 1
        elif self.DHA_angle <= 90: return (90 - self.DHA_angle)/30.0
        else: return 0

    def score(self):
        return self.dist_score()*self.angle_score()


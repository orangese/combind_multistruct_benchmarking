import math
from schrodinger.structutils.measure import measure_distance, measure_bond_angle

class HBond_Container:
    def __init__(self, lig, ind, settings):
        self.lig = lig
        self.ind = ind
        self.all_hdon = {}
        self.all_hacc = {}
        self.settings = settings

    def add_residue(self, resnum, res_st):
        self.all_hdon[resnum] = []
        self.all_hacc[resnum] = []

        for res_atom in res_st.hdon:
            for lig_atom in self.lig.hacc:
                for hydrogen in res_atom.bonded_atoms:
                    if hydrogen.element != 'H': continue
                    hbond = HBond(res_atom, lig_atom, hydrogen, resnum, True, self.settings)
                    if hbond.score():
                        self.all_hdon[resnum] += [hbond]

        for res_atom in res_st.hacc:
            for lig_atom in self.lig.hdon:
                for hydrogen in lig_atom.bonded_atoms:
                    if hydrogen.element != 'H': continue
                    hbond = HBond(lig_atom, res_atom, hydrogen, resnum, False, self.settings)
                    if hbond.score():
                        self.all_hacc[resnum] += [hbond]

    def filter_int(self):
        # enforces 1 hbond per h
        unique_h = {}
        for resnum in set(list(self.all_hdon.keys()) + list(self.all_hacc.keys())):
            for hb in self.all_hdon.get(resnum, []) + self.all_hacc.get(resnum, []):
                key = (hb.h.resnum, hb.h.index, hb.d.index)
                if key not in unique_h or unique_h[key].score() < hb.score():
                    unique_h[key] = hb

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
        for r, hb_list in self.all_hdon.items():
            for hb in hb_list:
                key = (self.ind[0], r, '')
                if key not in all_scores: all_scores[key] = 0
                all_scores[key] += hb.score()
        for r, hb_list in self.all_hacc.items():
            for hb in hb_list:
                key = (self.ind[1], r, '')
                if key not in all_scores: all_scores[key] = 0
                all_scores[key] += hb.score()
        return all_scores

class HBond:
    def __init__(self, donor, acceptor, h, r_ind, resIsHDonor, settings): # D - H ... A - X
        self.d = donor # h donor, \in {N,O}
        self.a = acceptor # h acceptor, \in {N,O}
        self.h = h # covalently bound to the donor
        self.resIsHDonor = resIsHDonor
        self.r_ind = r_ind
        self.settings = settings

        self.dist = measure_distance(h, acceptor) # angstroms
        self.DHA_angle = 180 - measure_bond_angle(self.d, self.h, self.a) # degrees

    def dist_score(self):
        if self.dist <= self.settings['hbond_dist_opt']:
            return 1
        elif self.dist <= self.settings['hbond_dist_cut']:
            return ((self.settings['hbond_dist_cut'] - self.dist)
                    / (self.settings['hbond_dist_cut'] - self.settings['hbond_dist_opt']))
        else:
            return 0

    def angle_score(self):
        if self.DHA_angle <= self.settings['hbond_angle_opt']:
            return 1
        elif self.DHA_angle <= self.settings['hbond_angle_cut']:
            return ((self.settings['hbond_angle_cut'] - self.DHA_angle)
                    / (self.settings['hbond_angle_cut'] - self.settings['hbond_angle_opt']))
        else:
            return 0

    def score(self):
        return self.dist_score() * self.angle_score()

    def __str__(self):
        return 'Donor: {}:{}, Acceptor: {}:{}, Score: {}'.format(self.d.resnum, self.d.pdbname,
                                                                 self.a.resnum, self.a.pdbname,
                                                                 self.score())


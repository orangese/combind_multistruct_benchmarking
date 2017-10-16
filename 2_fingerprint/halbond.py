from schrodinger.structutils.measure import measure_distance, measure_bond_angle
import math

class HalBond_Container:
    def __init__(self, lig_st):
        self.lig_st = lig_st
        self.all_haldon = {}
        self.all_halacc = {}

    def add_residue(self, resnum, res_st):
        self.all_haldon[resnum] = []
        self.all_halacc[resnum] = []
        for res_atom in res_st.atom:
            for lig_atom in self.lig_st.atom:
                if self.valid_donor(res_atom) and self.valid_acceptor(lig_atom):
                    self.all_haldon[resnum].extend([HalBond(res_atom, lig_atom, n, True) for n in res_atom.bonded_atoms if self.valid_hal(n)])
                if self.valid_donor(lig_atom) and self.valid_acceptor(res_atom):
                    self.all_halacc[resnum].extend([HalBond(lig_atom, res_atom, n, False) for n in lig_atom.bonded_atoms if self.valid_hal(n)])
    
    def valid_donor(self, atom):
        return atom.element in ('N', 'C') or self.valid_hal(atom)

    def valid_hal(self, atom):
        return atom.element in ('Cl', 'Br', 'I') and atom.formal_charge >= 0

    def valid_acceptor(self, atom): # h bond helper #2
        if atom.element not in ('N', 'O'):
            return False
        return atom.formal_charge <= 0

    def filter_int(self):
        # enforces 1 hbond per h
        #return
        res_h = {}
        for r in self.all_haldon:
            for hb in self.all_haldon[r]:
                unique_index = (r, hb.hal.index)
                if unique_index not in res_h or res_h[unique_index].score() < hb.score():
                    res_h[unique_index] = hb
        self.all_haldon = {r:[hb for (hb_r,hb_h), hb in res_h.items() if hb_r == r and hb.score() > 0.05] for r in [x[0] for x in res_h]}
        self.all_haldon = {r:self.all_haldon[r] for r in self.all_haldon if len(self.all_haldon[r]) > 0}

        lig_h = {}
        for r in self.all_halacc:
            for hb in self.all_halacc[r]:
                unique_index = (r, hb.hal.index)
                if unique_index not in lig_h or lig_h[unique_index].score() < hb.score():
                    lig_h[unique_index] = hb
        self.all_halacc = {r:[hb for (hb_r, hb_h), hb in lig_h.items() if hb_r == r and hb.score() > 0.05] for r in [x[0] for x in lig_h]}
        self.all_halacc = {r:self.all_halacc[r] for r in self.all_halacc if len(self.all_halacc[r]) > 0}        

    def all_res(self):
        return self.all_halacc.keys() + [r for r in self.all_haldon if r not in self.all_halacc]

    def score(self):
        return {
            r : [
                sum([hb.score() for hb in self.all_haldon.get(r, []) ]),
                sum([hb.score() for hb in self.all_halacc.get(r, []) ]),
                #sum([hb.new_score() for hb in self.all_haldon.get(r, []) ]),
                #sum([hb.new_score() for hb in self.all_halacc.get(r, []) ])
            ] for r in self.all_res()
        }

    def __str__(self):
        hdon = 'Hal Donors: \n' + ''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_haldon[r]]) for r in sorted(self.all_haldon.keys())])
        hacc = 'Hal Acceptors: \n' + ''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_halacc[r]]) for r in sorted(self.all_halacc.keys())])
        return hdon + '\n' + hacc

class HalBond:
    def __init__(self, donor, acceptor, hal, resIsHalDonor): # D - H ... A - X
        self.d = donor # h donor, \in {N,O}
        self.a = acceptor # h acceptor, \in {N,O}
        self.hal = hal # covalently bound to the donor
        self.resIsHalDonor = resIsHalDonor

        self.dist = measure_distance(hal, acceptor) # angstroms
        self.DHA_angle = 180 - measure_bond_angle(self.d, self.hal, self.a)#self.get_DHA_angle() # degrees

        self.charge_score = 1.2 if (self.d.formal_charge > 0 or self.a.formal_charge < 0) else 1

    def new_dist_score(self):

        radii = {'C':1.7, 'Cl':1.75, 'I':1.98, 'Br':1.85, 'O':1.52, 'N':1.55}
        min_r = radii[self.hal.element] + radii[self.a.element]

        if self.dist <= min_r*1.2: return 1
        elif self.dist <= min_r*1.5: return (1.5*min_r - self.dist)/(0.3*min_r)
        else: return 0

    def new_angle_score(self):
        if self.DHA_angle <= 45: return 1
        elif self.DHA_angle <= 90: return (90 - self.DHA_angle)/45.0
        else: return 0
    
    def score(self):
        return self.charge_score*self.new_dist_score()*self.new_angle_score()

    def __str__(self):
        d_str = '\n+donor: ' + str(self.d.index) + ' ' + self.d.element
        a_str = '\n+acceptor: ' + str(self.a.index) + ' ' + self.a.element
        h_str = '\n+halogen: ' + str(self.hal.index) + ' ' + self.hal.element
        dist  = '\n->dist: ' + str(self.dist)
        a1 = '\n->DHA: ' + str(self.DHA_angle)
        chrg = '\n->charge score: ' + str(self.charge_score) + '\n'
        return '\nHalBond score: ' + str(self.score()) + d_str + h_str + a_str + dist + a1 + chrg


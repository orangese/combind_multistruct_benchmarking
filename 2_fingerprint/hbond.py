from schrodinger.structutils.measure import measure_distance, measure_bond_angle
import math
class HBond_Container:
    def __init__(self, lig_st):
        self.lig_st = lig_st
        self.all_hdon = {}
        self.all_hacc = {}

    def add_residue(self, resnum, res_st):
        self.all_hdon[resnum] = []
        self.all_hacc[resnum] = []
        for res_atom in res_st.atom:
            for lig_atom in self.lig_st.atom:
                if self.valid_donor(res_atom) and self.valid_acceptor(lig_atom):
                    self.all_hdon[resnum].extend([HBond(res_atom, lig_atom, n, True) for n in res_atom.bonded_atoms if n.element == 'H'])
                if self.valid_donor(lig_atom) and self.valid_acceptor(res_atom):
                    self.all_hacc[resnum].extend([HBond(lig_atom, res_atom, n, False) for n in lig_atom.bonded_atoms if n.element == 'H'])
    
    def valid_donor(self, atom):
        return (atom.element in ('N', 'O')) and ('H' in [n.element for n in atom.bonded_atoms]) and atom.formal_charge >= 0

    def valid_acceptor(self, atom): # h bond helper #2
        if atom.element not in ('N', 'O'):
            return False
        return atom.formal_charge <= 0

    def filter_int(self):
        # enforces 1 hbond per h
        
        res_h = {}
        for r in self.all_hdon:
            for hb in self.all_hdon[r]:
                unique_index = (r, hb.h.index)
                if unique_index not in res_h or res_h[unique_index].score() < hb.score():
                    res_h[unique_index] = hb
        self.all_hdon = {r:[hb for (hb_r,hb_h), hb in res_h.items() if hb_r == r and hb.score() > 0.05] for r in [x[0] for x in res_h]}
        self.all_hdon = {r:self.all_hdon[r] for r in self.all_hdon if len(self.all_hdon[r]) > 0}

        lig_h = {}
        for r in self.all_hacc:
            for hb in self.all_hacc[r]:
                unique_index = (r, hb.h.index)
                if unique_index not in lig_h or lig_h[unique_index].score() < hb.score():
                    lig_h[unique_index] = hb
        self.all_hacc = {r:[hb for (hb_r, hb_h), hb in lig_h.items() if hb_r == r and hb.score() > 0.05] for r in [x[0] for x in lig_h]}
        self.all_hacc = {r:self.all_hacc[r] for r in self.all_hacc if len(self.all_hacc[r]) > 0}        

    def all_res(self):
        return self.all_hacc.keys() + [r for r in self.all_hdon if r not in self.all_hacc]

    def score(self):
        return {
            r : [
                #sum([hb.score() for hb in self.all_hdon.get(r, []) ]),
                #sum([hb.score() for hb in self.all_hacc.get(r, []) ]),
                sum([hb.new_score() for hb in self.all_hdon.get(r, []) ]),
                sum([hb.new_score() for hb in self.all_hacc.get(r, []) ])
            ] for r in self.all_res()
        }

    def __str__(self):
        hdon = 'H Donors: \n' + ''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_hdon[r]]) for r in sorted(self.all_hdon.keys())])
        hacc = 'H Acceptors: \n' + ''.join([str(r) + ''.join(['\n'+str(i) for i in self.all_hacc[r]]) for r in sorted(self.all_hacc.keys())])
        return hdon + '\n' + hacc

class HBond:
    def __init__(self, donor, acceptor, h, resIsHDonor): # D - H ... A - X
        self.d = donor # h donor, \in {N,O}
        self.a = acceptor # h acceptor, \in {N,O}
        self.h = h # covalently bound to the donor
        self.resIsHDonor = resIsHDonor

        self.dist = measure_distance(h, acceptor) # angstroms
        self.DHA_angle = 180 - measure_bond_angle(self.d, self.h, self.a)#self.get_DHA_angle() # degrees
        #self.HAX_angle = self.get_HAX_angle()

        self.charge_score = 1.2 if (self.d.formal_charge > 0 or self.a.formal_charge < 0) else 1

    #def get_DHA_angle(self):
    #    return measure_bond_angle(self.d, self.h, self.a)
        #return math.fabs(180 - func.angle_between_three_points(self.d.coordinates,
        #                            self.h.coordinates, self.a.coordinates) * 180 / math.pi)

    #def get_HAX_angle(self):
        #x = [atom.coordinates for atom in self.a.connected_atoms]
        #a = self.a.coordinates
        #unit_vectors = [func.return_normalized_vector(func.vector_subtraction(p, a)) for p in x]

        #if len(unit_vectors) == 3:
        #    unit_vectors = [
        #        func.return_normalized_vector(func.CrossProduct(unit_vectors[0], unit_vectors[1])),
        #        func.return_normalized_vector(func.CrossProduct(unit_vectors[1], unit_vectors[2])),
        #        func.return_normalized_vector(func.CrossProduct(unit_vectors[2], unit_vectors[0]))
        #    ]

        #    h_loc = func.vector_subtraction(self.h.coordinates, a)
        #    neg_v = func.vector_subtraction(Point(0,0,0), unit_vectors[0])
        #    if h_loc.dist_to(unit_vectors[0]) < h_loc.dist_to(neg_v):
        #        unit_vectors[0] = neg_v
        #    for i in range(1,3):
        #        (opt1, opt2) = (unit_vectors[i], func.vector_subtraction(Point(0,0,0), unit_vectors[i]))
        #        if unit_vectors[0].dist_to(opt2) < unit_vectors[0].dist_to(opt1):
        #            unit_vectors[i] = opt2

        #x_coords = func.vector_addition(func.vector_average(unit_vectors), a)

        #return math.fabs(180 - func.angle_between_three_points(self.h.coordinates, a, x_coords)*180/math.pi)

    def angle_score(self, isHAX):
        if isHAX and self.a.element == 'O': return 1
        #    return 1/(1+math.exp((self.HAX_angle-100)/10))
        elif isHAX: return 1
        #    return 1/(1+math.exp((self.HAX_angle-60)/10))
        return 1/(1+math.exp((self.DHA_angle-60)/10))

    def new_dist_score(self):
        if self.dist <= 2: return 1
        elif self.dist <= 3: return (3 - self.dist)
        else: return 0

    def new_angle_score(self):
        #if self.a.element == 'O':
        if self.DHA_angle <= 45: return 1
        elif self.DHA_angle <= 90: return (90 - self.DHA_angle)/45.0
        else: return 0

    def dist_score(self):
        return 1/(1+math.exp(3*(self.dist-2.6)))

    def score(self):
        return self.dist_score()*self.angle_score(False)*self.angle_score(True)*self.charge_score

    def new_score(self):
        return self.charge_score*self.new_dist_score()*self.new_angle_score()

    def __str__(self):
        d_str = '\n+donor: ' + str(self.d.index) + ' ' + self.d.element
        a_str = '\n+acceptor: ' + str(self.a.index) + ' ' + self.a.element
        h_str = '\n+hydrogen: ' + str(self.h.index)
        dist  = '\n->dist: ' + str(self.dist)
        a1 = '\n->DHA: ' + str(self.DHA_angle)
        #a2 = '\n->HAX: ' + str(self.HAX_angle)
        chrg = '\n->charge score: ' + str(self.charge_score) + '\n'
        return '\nHBond score: ' + str(self.score()) +','+ str(self.new_score()) + d_str + h_str + a_str + dist + a1 + chrg


import math
import math_functions as func
from point_atom import Point, Atom

def valid_donor(atom): # h bond helper #1
    return (atom.element in ('N', 'O')) and ('H' in [n.element for n in atom.connected_atoms])

def valid_acceptor(atom): # h bond helper #2
    if atom.element not in ('N', 'O'):
        return False
    return atom.formal_charge <= 0
    
def valid_sb(atom1, atom2):
    q1 = atom1.formal_charge
    q2 = atom2.formal_charge
    
    local_q1 = q1 + sum([n.formal_charge for n in atom1.connected_atoms])
    local_q2 = q2 + sum([n.formal_charge for n in atom2.connected_atoms])

    # (1) must have opposite formal charges
    # (2) formal charges must not be (approximately) cancelled by neighboring atoms
    return q1*q2 < 0 and local_q1*local_q2 < 0

class HBond:
    def __init__(self, donor, acceptor, h, resIsHDonor): # D - H ... A - X
        self.d = donor # h donor, \in {N,O}
        self.a = acceptor # h acceptor, \in {N,O}
        self.h = h # covalently bound to the donor
        self.resIsHDonor = resIsHDonor

        self.dist = h.dist_to(acceptor) # angstroms
        self.DHA_angle = self.get_DHA_angle() # degrees
        self.HAX_angle = self.get_HAX_angle()

    def get_DHA_angle(self):
        return math.fabs(180 - func.angle_between_three_points(self.d.coordinates, 
                                    self.h.coordinates, self.a.coordinates) * 180 / math.pi)

    def get_HAX_angle(self):
        x = [atom.coordinates for atom in self.a.connected_atoms]
        a = self.a.coordinates
        unit_vectors = [func.vector_addition(func.return_normalized_vector(func.vector_subtraction(p, a)), a) for p in x]
        
        if len(unit_vectors) == 3:
            unit_vectors = [
                func.return_normalized_vector(func.CrossProduct(unit_vectors[0], unit_vectors[1])),
                func.return_normalized_vector(func.CrossProduct(unit_vectors[1], unit_vectors[2])),
                func.return_normalized_vector(func.CrossProduct(unit_vectors[2], unit_vectors[0]))
            ]

            for (i, v) in enumerate(unit_vectors):
                neg_v = func.vector_subtraction(Point(0,0,0), v)
                if self.h.coordinates.dist_to(v) < self.h.coordinates.dist_to(neg_v):
                    unit_vectors[i] = neg_v

        x_coords = func.vector_average(unit_vectors)

        return math.fabs(180 - func.angle_between_three_points(self.h.coordinates, a, x_coords)*180/math.pi)

    def score(self):
        dist_score = 1/(1+math.exp(4*(self.dist-2.6)))
        DHA_score  = 1/(1+math.exp((self.DHA_angle-60)/10))
        HAX_score  = 1/(1+math.exp((self.HAX_angle-60)/10))
        return dist_score*DHA_score*HAX_score

    def __str__(self):
        d_str = '\n+donor: ' + str(self.d.atom_id) + ' ' + self.d.element
        a_str = '\n+acceptor: ' + str(self.a.atom_id) + ' ' + self.a.element
        h_str = '\n+hydrogen: ' + str(self.h.atom_id)
        dist  = '\n->dist: ' + str(self.dist)
        a1 = '\n->DHA: ' + str(self.DHA_angle)
        a2 = '\n->HAX: ' + str(self.HAX_angle)
        return '\nHBond score: ' + str(self.score()) + d_str + h_str + a_str + dist + a1 + a2

class LJ: # Lennard Jones potential
    def __init__(self):
        # values taken from Autodock 3.0.5 user guide
        # http://autodock.scripps.edu/faqs-help/manual/autodock-3-user-s-guide/AutoDock3.0.5_UserGuide.pdf
        # (C_{12} [kcal mol^-1 A^12], C_6 [kcal mol^-1 A^6])
        self.params = {
            'C':(2516582.400, 1228.800000),
            'N':(540675.281, 588.245000),
            'O':(230584.301, 429.496730),
            'S':(3355443.200, 1638.400000),
            'H':(81.920, 2.560000)
        }
        self.total_score = 0

    # add_score is called on all pairs of atoms (see residue.py)
    def add_score(self, atom1, atom2):
        (C12_1, C6_1) = self.params.get(atom1.element,(0,0))
        (C12_2, C6_2) = self.params.get(atom2.element,(0,0))
        r = atom1.dist_to(atom2)
        # count only favorable (negative, vdW) interactions 
        self.total_score += min(0, (C12_1*C12_2)**0.5/r**12 - (C6_1*C6_2)**0.5/r**6)

    def score(self): 
        return self.total_score

    def __str__(self):
        return '\nLJ potential: ' + str(self.score())

class SaltBridge:
    def __init__(self, atom1, atom2):
        # atom 1 and atom 2 must have opposite formal charges (see residue.py)
        self.atom1 = atom1
        self.atom2 = atom2
        self.r = atom1.dist_to(atom2)

    def score(self):
        # scales with 1/r (as in electric potential energy) and is capped at 2
        return min(2, 4.0/self.r)

    def __str__(self):
        at1 = '\n+atom1: ' + str(self.atom1.atom_id) + ' ' + self.atom1.element + '\n++charge: ' + str(self.atom1.formal_charge)
        at2 = '\n+atom2: ' + str(self.atom2.atom_id) + ' ' + self.atom2.element + '\n++charge: ' + str(self.atom2.formal_charge)
        return '\nSB score: ' + str(self.score()) + at1 + at2 + '\n+r: ' + str(self.r)


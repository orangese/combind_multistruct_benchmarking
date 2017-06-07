import math
import math_functions as func

class HBond:
    def __init__(self, donor, acceptor, h):
        self.donor = donor
        self.acceptor = acceptor
        self.h = h

    def atoms(self):
        return [self.donor, self.acceptor, self.h]

    def dist(self):
        return self.donor.dist_to(self.acceptor)

    def angle(self):
        return 180 / math.pi * func.angle_between_three_points(self.donor.coordinates,
                                               self.h.coordinates, self.acceptor.coordinates)
    def stat(self):
        return self.dist()

class Contact:
    def __init__(self, res_atom, lig_atom):
        self.res_atom = res_atom
        self.lig_atom = lig_atom

    def atoms(self):
        return [self.res_atom, self.lig_atom]

    def dist(self):
        return self.res_atom.dist_to(self.lig_atom)

    def coulumbic(self):
        return 100 * float((self.lig_atom.charge * self.res_atom.charge)) / (self.dist() ** 2)

    def hydrophobic(self):
        if self.lig_atom.charge < .1 and abs(self.res_atom.charge) < .2:
            return 1 / (self.dist() ** 2)
        return 0

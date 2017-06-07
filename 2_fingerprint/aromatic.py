from point_atom import Point
import math_functions as func
import math

class Aromatic:
    """
    An aromotic stores information about a given aromatic ring
    There is no checking done to see if a ring actually is aromatic
    This is assumed to be done correctly in Ligand or Residue
    """
    def __init__(self, atoms):
        self.atoms = atoms
        self.center = self._find_center()
        self.radius = self._find_radius()
        self.plane_coeff = self._find_plane_coeff()

    def _find_center(self):
        total = float(len(self.atoms))
        x_sum, y_sum, z_sum = 0.0, 0.0, 0.0
        for atom in self.atoms:
            x_sum += atom.coordinates.x
            y_sum += atom.coordinates.y
            z_sum += atom.coordinates.z
        return Point(x_sum / total, y_sum / total, z_sum / total)

    def _find_radius(self):
        return max(self.center.dist_to(atom.coordinates) for atom in self.atoms)

    def _find_plane_coeff(self):
        A = self.atoms[0].coordinates
        B = self.atoms[2].coordinates
        C = self.atoms[4].coordinates
        AB = func.vector_subtraction(B,A)
        AC = func.vector_subtraction(C,A)
        ABXAC = func.CrossProduct(AB,AC)
        # formula for plane will be ax + by + cz = d
        x1 = self.center.x
        y1 = self.center.y
        z1 = self.center.z
        a = ABXAC.x
        b = ABXAC.y
        c = ABXAC.z
        d = a*x1 + b*y1 + c*z1
        return (a, b, c, d)

TOLERANCE = 30
PIPI_DIST = 4.0
T_DIST = 4.4
OFFSET = 3
class Pi:
    def __init__(self, aromatic1, aromatic2):
        self.aromatic1 = aromatic1
        self.aromatic2 = aromatic2

    def atoms(self):
        return self.aromatic1.atoms + self.aromatic2.atoms

    def angle(self):
        aromatic1_norm_vector = Point(self.aromatic1.plane_coeff[0], self.aromatic1.plane_coeff[1], self.aromatic1.plane_coeff[2])
        aromatic2_norm_vector = Point(self.aromatic2.plane_coeff[0], self.aromatic2.plane_coeff[1], self.aromatic2.plane_coeff[2])
        return func.angle_between_points(aromatic1_norm_vector, aromatic2_norm_vector) * 180.0/math.pi

    def offset1(self):
        return func.project_point_onto_plane(self.aromatic1.center, self.aromatic2.plane_coeff).dist_to(self.aromatic2.center)

    def offset2(self):
        return func.project_point_onto_plane(self.aromatic2.center, self.aromatic1.plane_coeff).dist_to(self.aromatic1.center)

    #Minimum projected offset between aromatic rings
    def pi_pi_offset(self):
        return min(self.offset1(), self.offset2()) - self.aromatic1.radius - self.aromatic2.radius

    def pi_pi(self):
        if self.aromatic1.center.dist_to(self.aromatic2.center) > PIPI_DIST: return False
        if self.pi_pi_offset() < OFFSET and (abs(self.angle() - 180) < TOLERANCE or abs(self.angle()) < TOLERANCE):
            return self.pi_pi_offset(), self.dist(), self.angle()
        return False

    def pi_t(self):
        if self.aromatic1.center.dist_to(self.aromatic2.center) > T_DIST: return False
        if abs(self.angle() - 90) > TOLERANCE and abs(self.angle() - 270) > TOLERANCE: return False
        if self.offset1() < self.aromatic2.radius + OFFSET:
            return self.dist()
        elif self.offset2() < self.aromatic1.radius + OFFSET:
            return self.dist()
        return False

    def dist(self):
        return self.aromatic1.center.dist_to(self.aromatic2.center)

class PiCation:
    def __init__(self, aro, charged):
        self.aro = aro
        self.charged = charged

    def atoms(self):
        return self.aro.atoms + [self.charged]

    def dist(self):
        return self.aro.center.dist_to(self.charged.coordinates)

    def stat(self):
        return self.dist()

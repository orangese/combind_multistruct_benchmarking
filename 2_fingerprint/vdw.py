import math
from schrodinger.structutils.measure import measure_distance

class VDW_Container: # Lennard Jones potential
    def __init__(self, lig_st):
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
        self.total = {}
        self.total_other = {}
        self.lig_st = lig_st

    def add_residue(self, resnum, res_st):
        for res_atom in res_st.atom:
            for lig_atom in self.lig_st.atom:
                self.score_atom_pair(res_atom, lig_atom, resnum)
    def filter_int(self): pass
    def score_atom_pair(self, atom1, atom2, resnum):
        (C12_1, C6_1) = self.params.get(atom1.element,self.params['C'])
        (C12_2, C6_2) = self.params.get(atom2.element,self.params['C'])
        (C12, C6) = ((C12_1*C12_2)**0.5, (C6_1*C6_2)**0.5)

        r = measure_distance(atom1, atom2)
        
        if resnum not in self.total:
            self.total[resnum] = 0
            self.total_other[resnum] = 0

        # count only favorable (negative, vdW) interactions 
        self.total[resnum] += C12/r**12 - C6/r**6

        r_m = (2*C12/C6)**(0.16666667)
        eps = -(C6**2)/(4*C12)
        if abs(r-r_m) <= 0.5:
            self.total_other[resnum] += eps
        elif r > (r_m + 0.5):
            self.total_other[resnum] += C12/(r-0.5)**12 - C6/(r-0.5)**6
        elif r < (r_m - 0.5):
            self.total_other[resnum] += C12/(r+0.5)**12 - C6/(r+0.5)**6
    
    def score(self): 
        return {
            r : [
                max(0, -1*self.total[r]),
                max(0, -1*self.total_other[r])
            ] for r in self.total
        }

    #def __str__(self):
    #    return '\nLJ potential: ' + str(self.score()[r])

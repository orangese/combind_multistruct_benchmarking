from schrodinger.structutils.measure import measure_distance

class Hydrophobic_Container:
    def __init__(self, lig_st):
        self.lig_st = lig_st

        self.radii = {'H': 1.2, 'C':1.7, 'F':1.47, 'Cl':1.75, 'I':1.98, 'Br':1.85}
        self.all_scores = {}

    def is_hydrophobic(self, atom):
        if atom.element == 'H' and [n.element for n in atom.bonded_atoms][0] in self.radii:
            return True
        return atom.element in self.radii

    def add_residue(self, resnum, res_st):
        self.all_scores[resnum] = 0
        for res_atom in res_st.atom:
            if not self.is_hydrophobic(res_atom): continue
            max_coverage = 0
            for lig_atom in self.lig_st.atom:
                if not self.is_hydrophobic(lig_atom): continue
                vdw_radius = self.radii[res_atom.element] + self.radii[lig_atom.element]
                r = measure_distance(res_atom, lig_atom)
                if r <= 1.25*vdw_radius:
                    max_coverage += 1
                elif r <= 1.75*vdw_radius:
                    max_coverage += (1.75*vdw_radius - r)/(0.5*vdw_radius)
            self.all_scores[resnum] += max_coverage

    def filter_int(self):
        pass

    def score(self):
        return { r : [ self.all_scores[r] ] for r in self.all_scores }

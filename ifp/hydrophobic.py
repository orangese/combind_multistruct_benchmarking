from schrodinger.structutils.measure import measure_distance

class Hydrophobic_Container:
    def __init__(self, lig, ind, settings):
        self.radii = {'H': 1.2, 'C':1.7, 'F':1.47, 'Cl':1.75, 'I':1.98, 'Br':1.85}

        self.lig = lig
        self.ind = ind
        self.settings = settings

        self.all_scores = {}

    def add_residue(self, resnum, res):
        for res_atom in res.nonpolar1 | res.nonpolar2:
            for lig_atom in self.lig.nonpolar1 | self.lig.nonpolar2:

                vdw_radius = self.radii[res_atom.element] + self.radii[lig_atom.element]
                r = measure_distance(res_atom, lig_atom)
                if r >= self.settings['contact_scale_cut']*vdw_radius: continue
   
                score = 0
                if r <= self.settings['contact_scale_opt']*vdw_radius:
                    score = 1
                else:
                    score = ((self.settings['contact_scale_cut']*vdw_radius - r)
                             / ((self.settings['contact_scale_cut'] - self.settings['contact_scale_opt'])
                                * vdw_radius))
            
                if lig_atom in self.lig.nonpolar1 and res_atom in res.nonpolar1:
                    key1 = (self.ind[0], resnum, '')
                    self.all_scores[key1] = self.all_scores.get(key1, 0) + score
                if lig_atom in self.lig.nonpolar2 and res_atom in res.nonpolar2:
                    key2 = (self.ind[1], resnum, '')
                    self.all_scores[key2] = self.all_scores.get(key2, 0) + score
    
    def filter_int(self):
        pass

    def score(self):
        return self.all_scores

    def raw(self):
        return {}

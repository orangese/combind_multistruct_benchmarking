from schrodinger.structutils.measure import measure_distance

class SB_Container:
    def __init__(self, lig, ind, settings):
        self.lig = lig
        self.ind = ind
        self.settings = settings
        self.all_sb = {}

    def add_residue(self, resnum, res):
        self.all_sb[resnum] = []
        for res_atom in res.chrg:
            for lig_atom in self.lig.chrg:
                sb = SB(res_atom, lig_atom, self.settings)
                if sb.is_valid():
                    self.all_sb[resnum].append(sb)

    def filter_int(self):
        pass

    def score(self):
        all_scores = {}
        for r, sb_list in self.all_sb.items():
            for sb in sb_list:
                k = (self.ind[0], r,'')
                if k not in all_scores: all_scores[k] = 0
                all_scores[k] += sb.score()
        return all_scores

class SB:
    def __init__(self, res_atom, lig_atom, settings):
        self.res_atom = res_atom
        self.lig_atom = lig_atom
        self.settings = settings
        self.dist = None

    def is_valid(self):
        """
        Checks if atoms have opposite formal charges and are within
        range to form a saltbridge.
        """
        return ((self.res_atom.formal_charge*self.lig_atom.formal_charge) < 0
                and self.score())

    def score(self):
        if self.dist is None:
            self.dist = measure_distance(self.res_atom, self.lig_atom)
        
        if self.dist <= self.settings['sb_dist_opt']:
            return 1
        elif self.dist <= self.settings['sb_dist_cut']:
            return ((self.settings['sb_dist_cut'] - self.dist)
                    / (self.settings['sb_dist_cut'] - self.settings['sb_dist_opt']))
        else:
            return 0

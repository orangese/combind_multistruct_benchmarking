from schrodinger.structutils.measure import measure_distance
from shared_paths import shared_paths
params = shared_paths['ifp']

class SB:
    def __init__(self, res_atom, lig_atom):
        self.res_atom = res_atom
        self.lig_atom = lig_atom
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
        
        if self.dist <= params['sb_dist_opt']:
            return 1
        elif self.dist <= params['sb_dist_cut']:
            return ((params['sb_dist_cut'] - self.dist)
                    / (params['sb_dist_cut'] - params['sb_dist_opt']))
        else:
            return 0

class SB_Container:
    def __init__(self, lig, indices):
        self.lig = lig
        self.indices = indices
        self.all_sb = {}

    def add_residue(self, resnum, res):
        self.all_sb[resnum] = []
        for res_atom in res.chrg:
            for lig_atom in self.lig.chrg:
                sb = SB(res_atom, lig_atom)
                if sb.is_valid():
                    self.all_sb[resnum].append(sb)

    def filter_int(self):
        """
        Previously used to enforce one saltbridge per ligand formal charge,
        but on second thought, there doesn't seem to be a need to do this.
        """
        pass

    def score(self):
        """
        Computes final salt bridge fingerprints.
        """
        all_scores = {}
        for r, sb_list in self.all_sb.items():
            for sb in sb_list:
                sc = [sb.score()]
                for i, j in enumerate(self.indices):
                    k = (j, r,'')
                    if k not in all_scores: all_scores[k] = 0
                    all_scores[k] += sc[i]
        return all_scores

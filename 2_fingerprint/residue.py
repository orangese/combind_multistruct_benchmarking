from pdb import PDB
from interactions import HBond, LJ, SaltBridge, valid_donor, valid_acceptor, valid_sb

class Residue(PDB):
    """ A Residue is one of the standard AAs
    Its atoms should be named in the standard way.
    This class should be used to represent a subset of a larger pdb
    """
    def __init__(self, resname):
        self.resname = resname.strip()
        self.num = -1
        self.atoms = []
        self.aromatics = []

        self.score_thresh = 0.1

        self.all_interactions = []
        self.fp = [0,0,0,0]

    def copy_of(self):
        """
        Return fresh copy with no interaction information
        """
        new = Residue(self.resname)
        for atom in self.atoms:
            new.add_atom(atom.copy_of())
        new.assign()
        return new

    def add_atom(self, atom):
        if self.num != -1: assert self.num == atom.residue_id
        self.num = atom.residue_id
        self.atoms.append(atom)

    def all_atoms(self):
        return [atom for atom in self.atoms]

    def fingerprint(self, ligand):
        if sum(self.fp) != 0: return self.fp
        
        hbonds = {} # maps h id to hbond (one h bond per h)
        sbs = []
        lj = LJ()

        for res_atom in self.atoms:
            for lig_atom in ligand.all_atoms():

                hb_options = []
                if valid_donor(res_atom) and valid_acceptor(lig_atom):
                    hb_options.extend([HBond(res_atom, lig_atom, n, True) for n in res_atom.connected_atoms if n.element == 'H'])
                if valid_donor(lig_atom) and valid_acceptor(res_atom):
                    hb_options.extend([HBond(lig_atom, res_atom, n, False) for n in lig_atom.connected_atoms if n.element == 'H'])

                for hb in hb_options:
                    if hb.h.atom_id not in hbonds or hbonds[hb.h.atom_id].score() < hb.score():
                            hbonds[hb.h.atom_id] = hb
                
                if valid_sb(res_atom, lig_atom):
                    sbs.append(SaltBridge(res_atom, lig_atom))

                lj.add_score(res_atom, lig_atom)

        hbonds = [hbonds[h] for h in hbonds.keys() if hbonds[h].score() >= self.score_thresh]
        sbs = [i for i in sbs if i.score() >= self.score_thresh]

        self.fp[0] = sum([hb.score() for hb in hbonds if hb.resIsHDonor])
        self.fp[1] = sum([hb.score() for hb in hbonds if not hb.resIsHDonor])
        self.fp[2] = sum([sb.score() for sb in sbs])
        self.fp[3] = lj.score() if lj.score() >= self.score_thresh else 0

        self.all_interactions = hbonds + sbs
        
        return self.fp

    def assign(self):
        self._assign_bonds()
        self._assign_aromatic_rings()

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

    def get_interactions(self, ligand):
        hbonds = []
        sbs = []
        lj = LJ()

        for res_atom in self.atoms:
            for lig_atom in ligand.all_atoms():

                if valid_donor(res_atom) and valid_acceptor(lig_atom):
                    hbonds.extend([HBond(res_atom, lig_atom, n, True) for n in res_atom.connected_atoms if n.element == 'H'])
                if valid_donor(lig_atom) and valid_acceptor(res_atom):
                    hbonds.extend([HBond(lig_atom, res_atom, n, False) for n in lig_atom.connected_atoms if n.element == 'H'])

                if valid_sb(res_atom, lig_atom):
                    sbs.append(SaltBridge(res_atom, lig_atom))

                lj.add_score(res_atom, lig_atom)

        return hbonds, sbs, lj

    def assign(self):
        self._assign_bonds()
        self._assign_aromatic_rings()

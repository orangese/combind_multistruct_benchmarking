from point_atom import Atom
from pdb import PDB

class Ligand(PDB):
    """
    A ligand is a collection of atoms that are fully connecting
    They should represent a complete molecule
    AA residues do not fit this description and should be included as a Residue

    A ligand can represent part of a larger PDB as part of Receptor or be a full entry by itself

    A ligand does not have any further divisions
    """
    def __init__(self):
        self.atoms = []
        self.aromatics = []
        self.name = ''
        self.st = None

    def assign(self):
        self._assign_bonds()
        self._assign_aromatic_rings()

    def all_atoms(self):
        return self.atoms

    def debug_aromatics(self):
        return self.atoms + map(lambda x: x.center, self.aromatics)

    def add_atom(self, atom):
        self.name = atom.residue_id
        self.atoms.append(atom)

    def load_mae(self, st):
        self.st = st #Save a copy of the st in the class
        
        for atom in st.atom:
            cur_atom = Atom()
            cur_atom.from_schrod(atom)
            self.add_atom(cur_atom)
        self.assign()

    def molecular_weight(self):
        weight = 0
        for atom in self.atoms:
            a_w = self._atom_weight(atom.element)
            if a_w == 0: # atom type not found
                return atom.element
            weight += a_w
        return weight

    def diameter(self):
        maximum = 0
        for i1 in range(len(self.atoms)):
            for i2 in range(i1+1,len(self.atoms)):
                maximum = max(maximum, self.atoms[i1].dist_to(self.atoms[i2]))
        return maximum
                

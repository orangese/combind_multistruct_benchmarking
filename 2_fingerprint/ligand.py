from point_atom import Point, Atom
from aromatic import Aromatic
import math_functions as func
from pdb import PDB
import textwrap
import math

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

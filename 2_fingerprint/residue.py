from pdb import PDB
from aromatic import Aromatic, Pi, PiCation
from non_aro_interactions import HBond, Contact
import math_functions as func
from point_atom import Point
import math
from schrodinger.structutils.minimize import Minimizer
from schrodinger.structure import Structure

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

	self.h_donor = 0
	self.h_acceptor = 0
	self.electrostatics = 0
	self.lj = 0

        self.to_print = []

    def copy_of(self):
        """
        Return fresh copy with no interaction information
        """
        new = Residue(self.resname)
        for atom in self.atoms:
            new.add_atom(atom.copy_of())
        new.assign()
        return new

    def add_electrostatic_potential(self, score):
        self.electrostatics += score

    # not sure what to do about the repulsive part -- 0 for now
    def add_lj_potential(self, score):
        if score < 0: self.lj += score

    def debug_h(self, lig_a, res_a, h_a, is_donor, score):
        if score > 0.3:
	    self.to_print.append((lig_a.atom_id, res_a.atom_id, h_a.atom_id, is_donor, score))

    def add_h_bond(self, score, is_donor):
        if is_donor:
            self.h_donor += score
        else:
            self.h_acceptor += score

    def add_atom(self, atom):
        if self.num != -1: assert self.num == atom.residue_id
        self.num = atom.residue_id
        self.atoms.append(atom)

    def all_atoms(self):
        return [atom for atom in self.atoms]

    def interacting_atoms(self):
        return set([atom for interaction in self.interactions.values() for x in interaction for atom in x.atoms()])
    
    ## hbonds get scaled x10
    ## potentials get capped to +- 10
    def fingerprint(self):#, combinedStruct):
        return [
                   10*self.h_donor, 
                   10*self.h_acceptor, 
                   min(max(self.electrostatics, -10), 10), 
                   min(max(self.lj, -10), 10)
               ] # + self.to_print

    def debug_fp(self):
        return self.debug_h

    def assign(self):
        self._assign_bonds()
        self._assign_aromatic_rings()

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

        self.interactions = {
            'h_donor'     : [],
            'h_acceptor'  : [],
            'contact'     : [],
            'contact_orig': [],
            'pi_stack'    : [],
            'pi_t'        : [],
            'pi_cation'   : []
            }

    def copy_of(self):
        """
        Return fresh copy with no interaction information
        """
        new = Residue(self.resname)
        for atom in self.atoms:
            new.add_atom(atom.copy_of())
        new.assign()
        return new

    def add_h_donor(self, donor, acceptor):
        self.interactions['h_donor'] += [(donor, acceptor)]

    def add_h_acceptor(self, donor, acceptor):
        self.interactions['h_acceptor'] += [(donor, acceptor)]

    def add_pi_stack(self, atom_list1, atom_list2):
        self.interactions['pi_stack'] += [(atom_list1, atom_list2)]
        
    def add_pi_t(self, atom_list1, atom_list2):
        self.interactions['pi_t'] += [(atom_list1, atom_list2)]

    def add_contact(self, res_atom, lig_atom):
        self.interactions['contact'] += [(res_atom, lig_atom)]
        
    def add_contact_orig(self, res_atom, lig_atom):
        self.interactions['contact_orig'] += [Contact(res_atom, lig_atom)]
    
    def add_pi_cation(self, atom_list1, atom_list2):
        self.interactions['pi_cation'] += [(atom_list1, atom_list2)]

    def add_atom(self, atom):
        if self.num != -1: assert self.num == atom.residue_id
        self.num = atom.residue_id
        self.atoms.append(atom)

    def all_atoms(self):
        return [atom for atom in self.atoms]

    def interacting_atoms(self):
        return set([atom for interaction in self.interactions.values() for x in interaction for atom in x.atoms()])
    
    def fingerprint(self, combinedStruct):
        
        minimizer = Minimizer(struct=combinedStruct, max_steps = 0)
        
        #NOTE: getNBEnergyFromTwoAtomLists is zero-indexed!
        #NOTE: THE API IS RIDICULOUSLY DUMB - getNBEnergyFromTwoAtomLists returns (L-J, Columbic) NOT THE OTHER WAY AROUND
        #EVEN THOUGH THAT'S WHAT THE DOCUMENTATION SAID
        
        h_donor_energy = [minimizer.getNBEnergyForTwoAtomLists([one-1], [two-1])[1]
                          for one, two in self.interactions['h_donor']]
        #print("H DONOR:")
       	#print(self.interactions['h_donor'])
        #print(h_donor_energy)
        
        h_acceptor_energy = [minimizer.getNBEnergyForTwoAtomLists([one-1], [two-1])[1]
                             for one, two in self.interactions['h_acceptor']]
        #print("H ACCEPTOR:")
        #print(self.interactions['h_acceptor'])
        #print(h_acceptor_energy)
        
        pi_stack_energy = [minimizer.getNBEnergyForTwoAtomLists([x-1 for x in one], [x-1 for x in two])[1]
                           for one, two in self.interactions['pi_stack']]
        #print("PI STACK:")
        #print(self.interactions['pi_stack'])
        #print(pi_stack_energy)
        
        pi_t_energy = [minimizer.getNBEnergyForTwoAtomLists([x-1 for x in one], [x-1 for x in two])[1]
                       for one, two in self.interactions['pi_t']]
        #print("PI T:")
        #print(self.interactions['pi_t'])
        #print(pi_t_energy)
        
        pi_cation_energy = [minimizer.getNBEnergyForTwoAtomLists([x-1 for x in one], [x-1 for x in two])[1]
                            for one, two in self.interactions['pi_cation']]
        #print("PI CATION:")
        #print(self.interactions['pi_cation'])
        #print(pi_cation_energy)
        
        contact_energy = [minimizer.getNBEnergyForTwoAtomLists([one-1], [two-1])[1]
                          for one, two in self.interactions['contact']]
        #print("CONTACT:")
        #print(self.interactions['contact'])
        #print(contact_energy)

        #DO NOT REMOVE
        #So I checked and the atom mapping from original contacts -> new contacts is right, it's just that the 
        #getNBEnergyForTwoAtomLists function gives me a different answer
        #NOTE: The signs returned by contact orig and contact are the same -> good sanity check, just in diff units?
        #print("CONTACT ORIG:")
        #print([(y.atoms()[0].atom_id, y.atoms()[1].atom_id) for y in self.interactions['contact_orig']])
        #print([(y.atoms()[0].charge, y.atoms()[1].charge) for y in self.interactions['contact_orig']])
        #print(map(lambda x: x.coulumbic(), self.interactions['contact_orig']))

        return [sum(h_donor_energy), sum(h_acceptor_energy), sum(pi_stack_energy), sum(pi_t_energy), sum(pi_cation_energy),
               sum(contact_energy)]

    def assign(self):
        self._assign_bonds()
        self._assign_aromatic_rings()

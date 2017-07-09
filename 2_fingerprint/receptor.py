# BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob
# Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite 
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.

# The below is implemented by Joe Paggi and largely based on BINANA

from point_atom import Point, Atom
from pdb import PDB
from residue import Residue
from ligand import Ligand
import copy

class Receptor(PDB):
    def __init__(self):
        self.residues = {}
        self.ligands = {}
        self.st = None #This contains both the ligand and protein
        
        #For use in close_pdb only
        self.receptorOnlySt = None #This is a structure that contains only the protein
        self.atomIDToRes = {} #Calculated in interactions.py, based off the combined atom numbering
        self.receptorAtoms = [] #Indexes of receptor atoms in the combined structure
        self.ligandAtoms = [] #Indexes of ligand atoms in the combined structure
        self.origAtomToCombined = {} #We convert from ligand/protein specific numbering using str(Coordinate)->atom number

    def seperate_atoms(self, receptor, ligand):      
        recepAtoms = []
        ligAtoms = []
        recepCoords = [str(b.coordinates) for a in receptor.residues.values() for b in a.atoms]
        ligCoords = [str(a.coordinates) for a in ligand.atoms]
        
        for atom in self.st.atom:
            #Dict to convert from points -> numbering
            tempPoint = str(Point(atom.x, atom.y, atom.z))
            self.origAtomToCombined[tempPoint] = atom.index
            
            #Make a list of the receptor and ligand atoms
            if tempPoint in recepCoords:
                recepAtoms.append(atom.index)
            if tempPoint in ligCoords:
                ligAtoms.append(atom.index)
        
        self.receptorAtoms = recepAtoms
        self.ligandAtoms = ligAtoms
        self.receptorOnlySt = copy.copy(self.st)
        self.receptorOnlySt.deleteAtoms(ligAtoms) #Remove the ligand from the st, this is needed for some interaction methods
                
    def all_atoms(self):
        return [i for pdb in self.residues.values() + self.ligands.values() for i in pdb.all_atoms()]

    def fingerprint(self, ligand):
        """
        Let's deal with processing this in another class?? 
        Alternatively we could make a fp class and add a method 
        that lets us add things as we go???
        """
        return {num: resi.fingerprint(ligand) for num, resi in self.residues.items() if any(resi.fingerprint(ligand))}

    def assign(self):
        for residue in self.residues.values(): residue.assign()
        for ligand in self.ligands.values(): ligand.assign()

    def add_residue(self, residue):
        assert residue.num not in self.residues
        self.residues[residue.num] = residue

    def load_mae(self, st):
        self.st = st #Save a copy of the structure in the class
        
        for atom in st.atom:
            cur_atom = Atom()
            cur_atom.from_schrod(atom)
            if cur_atom.residue_name in self.protein_resnames:
                if cur_atom.residue_id not in self.residues: self.residues[cur_atom.residue_id] = Residue(cur_atom.residue_name)
                self.residues[cur_atom.residue_id].add_atom(cur_atom)
                self.atomIDToRes[cur_atom.atom_id] = cur_atom.residue_id
            else:
                if cur_atom.residue_id not in self.ligands: self.ligands[cur_atom.residue_id] = Ligand()
                self.ligands[cur_atom.residue_id].add_atom(cur_atom)
                self.atomIDToRes[cur_atom.atom_id] = cur_atom.residue_id
        self.assign()

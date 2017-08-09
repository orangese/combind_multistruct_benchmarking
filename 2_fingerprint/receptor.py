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
        used_h = {} # enforces 1 hbond per h
        used_fc = {} # 1 sb per ligand formal charge
        lj = {}
        for num, resi in self.residues.items():
            h, s, l = resi.get_interactions(ligand)
            for hb in h:
                if hb.h.atom_id not in used_h or used_h[hb.h.atom_id][0].score() < hb.score():
                    used_h[hb.h.atom_id] = (hb, num)
            for sb in s:
                if sb.lig_atom.atom_id not in used_fc or used_fc[sb.lig_atom.atom_id][0].score() < sb.score():
                    used_fc[sb.lig_atom.atom_id] = (sb, num)
            lj[num] = l

        hbs = {num:[] for num in self.residues}
        sbs = {num:[] for num in self.residues}
        score_thresh = 0.05
        for h in used_h:
            hb, r = used_h[h]
            if hb.score() >= score_thresh: 
                hbs[r].append(hb)
        for fc in used_fc:
            sb, r = used_fc[fc]
            if sb.score() >= score_thresh:
                sbs[r].append(sb)            
            
        #score_thresh = 0.05
        fp = {r:[0,0,0,0,0] for r in self.residues}
        for r in self.residues:
            fp[r][0] = sum([hb.score() for hb in hbs[r] if hb.resIsHDonor])# and hb.score() >= score_thresh])
            fp[r][1] = sum([hb.score() for hb in hbs[r] if not hb.resIsHDonor])# and hb.score() >= score_thresh])
            fp[r][2] = sum([sb.score() for sb in sbs[r]])# if abs(sb.score()) >= score_thresh])
            fp[r][3] = abs(lj[r].score()) if abs(lj[r].score()) >= score_thresh else 0
            fp[r][4] = abs(lj[r].other_score()) if abs(lj[r].other_score()) >= score_thresh else 0
            #fp.extend(fp_r)    
        #return {num: resi.fingerprint(ligand) for num, resi in self.residues.items() if any(resi.fingerprint(ligand))}
        self.int_per_res = {r: hbs[r] + sbs[r] for r in self.residues}
        return {r:fp[r] for r in fp if any(fp[r])}

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

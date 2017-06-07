#print BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob
# Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite 
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.

# The below is implemented by Joe Paggi and largely based on BINANA

from pdb import PDB
import math_functions as func
from point_atom import Point, Atom
from receptor import Receptor
from ligand import Ligand
import textwrap
import math
from schrodinger.structutils import interactions
from schrodinger.structutils.interactions import SaltBridgeFinder
from schrodinger.structutils.interactions.hbond import get_hydrogen_bonds
from schrodinger.structure import Structure
from schrodinger.structure import StructureReader
from schrodinger.structure import _StructureAtom
from schrodinger.structure import get_pyatom_from_cppatom

RESIDUE_THRESH = 5 #Value is in angstroms

class Interactions:
    """
    Class used to create a fingerprint given a receptor, ligand and parameters

    This class should not mutate the given receptor!!!!!!!!
    No issues mutating the ligand however.
    """
    def __init__(self, receptor, ligand, parameters):
        self.parameters = parameters
        self.ligand = ligand
        self.receptor = receptor
        
        #Receptor residues that are within RESIDUE_THRESH are added to a new Receptor object
        self.close_pdb = self.get_close_pdb(receptor.st, ligand.st)
        
        #Initialize the interactions that we care about in our fingerprint
        self.get_hbond()
        self.get_pi()
        self.get_pi_cation()
        self.get_contacts()
        
        #Calculate our fingerprint - definitions in Residue class
        self.fp = self.close_pdb.fingerprint()

    def get_fp(self):
        return self.fp

    def get_close_pdb(self, receptorStruct, ligandStruct):
        close_pdb = Receptor()
        
        #Only add residues that are within RESIDUE_THRESH of the ligand center coordinate
        for residue in self.receptor.residues.values():
            added = False
            for ligand_atom in self.ligand.atoms:
                for receptor_atom in residue.atoms:
                    if ligand_atom.coordinates.dist_to(receptor_atom.coordinates) < RESIDUE_THRESH:
                        assert not added
                        close_pdb.add_residue(residue.copy_of())
                        added = True
                        break
                if added: break
        close_pdb.ligands = self.receptor.ligands
        close_pdb.assign()
        
        #Calculate Atom ID -> Res, based on the merged structure
        close_pdb.st = self.receptor.st
        atomIDToRes = {}
        
        for atom in close_pdb.st.atom:
            cur_atom = Atom()
            cur_atom.from_schrod(atom)
            atomIDToRes[cur_atom.atom_id] = cur_atom.residue_id
        
        close_pdb.atomIDToRes = atomIDToRes
        
        close_pdb.seperate_atoms(self.receptor, self.ligand)
        return close_pdb

    '''
    Adds contacts between atoms that are < 4.5 A from each other
    '''
    def get_contacts(self):
        for lig_atom in self.ligand.all_atoms():
            for residue in self.close_pdb.residues.values():
                for res_atom in residue.atoms:
                    if res_atom.dist_to(lig_atom) > self.parameters['hydrophobic_dist_cutoff']: continue
                    residue.add_contact_orig(res_atom, lig_atom)
                    residue.add_contact(self.close_pdb.origAtomToCombined[str(res_atom.coordinates)],
                                        self.close_pdb.origAtomToCombined[str(lig_atom.coordinates)])

    '''
    Adds H donors/acceptors to the close_pdb object
    
    H donors/acceptors have to be bonded to [N, O] atoms, within dist_cutoff and within angle_cutoff
    
    If the hydrogen is in the receptor, then we say it's a donor
    If the hydrogen is in the ligand, then we say it's a acceptor
    
    Dist < hydrogen_bond_dist_cutoff
    Angle < hydrogen_bond_angle_cutoff
    '''
    def get_hbond(self):
        #The Schrodinger API lets us find H bonds using the combined structure numbering, no old->new ID conversion is needed
        hbonds = get_hydrogen_bonds(self.close_pdb.st, atoms1=self.close_pdb.receptorAtoms,
                                    atoms2=self.close_pdb.ligandAtoms)
        
        for atom1, atom2 in hbonds:
            atom1_num = atom1.index
            atom2_num = atom2.index
            
            resNum1 = self.close_pdb.atomIDToRes[atom1_num]
            resNum2 = self.close_pdb.atomIDToRes[atom2_num]
            
            if resNum1 in self.close_pdb.residues:
                self.close_pdb.residues[resNum1].add_h_donor(atom1_num, atom2_num)
                
            if resNum2 in self.close_pdb.residues:
                self.close_pdb.residues[resNum2].add_h_acceptor(atom1_num, atom2_num)
                
        '''
        for lig_atom in self.ligand.all_atoms():
            if lig_atom.element not in ('N', 'O'): continue

            for residue in self.close_pdb.residues.values():
                for res_atom in residue.atoms:
                    if res_atom.element not in ('N', 'O'): continue
                    dist = res_atom.dist_to(lig_atom)
                    if dist > self.parameters['hydrogen_bond_dist_cutoff']: continue
                    receptor_hydrogens = [atom for atom in list(res_atom.connected_atoms) if atom.element == 'H']
                    ligand_hydrogens = [atom for atom in list(lig_atom.connected_atoms) if atom.element == 'H']

                    for h in receptor_hydrogens:
                        angle = math.fabs(180 - func.angle_between_three_points(lig_atom.coordinates,
                                                                                h.coordinates, res_atom.coordinates) * 180 / math.pi)
                        if angle <= self.parameters['hydrogen_bond_angle_cutoff']:
                            residue.add_h_donor(res_atom, lig_atom, h)
                            
                    for h in ligand_hydrogens:
                        angle = math.fabs(180 - func.angle_between_three_points(lig_atom.coordinates,
                                                                                h.coordinates, res_atom.coordinates) * 180 / math.pi)
                        if angle <= self.parameters['hydrogen_bond_angle_cutoff']:
                            residue.add_h_acceptor(lig_atom, res_atom, h)
        '''

    '''
    Adds potential pi-pi stacks to the close_pdb object
    Potential pi-pi stacks are simply evaluated on center distances < pi_pi_interacting_dist_cutoff
    
    Note that validity of pi-pi stacks are evaluated in the residue class and are seperated into
    sandwich/parallel and T-shaped stacks in the feature vector
    '''
    def get_pi(self):
        #Schrodinger API needs seperate structures, do the old->new conversion
        #NOTE: THIS CONTAINS CUSTOM CODE THAT I (THOMAS) PUT IN TO RETURN ATOMS INSTEAD OF STRUCTS
        pipis = interactions.find_pi_pi_interactions(self.close_pdb.receptorOnlySt, struct2=self.ligand.st)
        
        for ringAtoms1, ringAtoms2, face_to_face in pipis:
            #NOTE: STRUCTURES RETURNED ARE STRUCTURE ATOMS: 
            #http://content.schrodinger.com/Docs/r2017-1/python_api/api/schrodinger.structure._StructureAtom-class.html
            
            struct1_atoms = []
            for x in ringAtoms1:
                x = get_pyatom_from_cppatom(x)
                #Convert from receptor and ligand specific numbering to the combined receptor + ligand numbering
                atomInd = self.close_pdb.origAtomToCombined[str(Point(_StructureAtom._getX(x), _StructureAtom._getY(x),
                                                                      _StructureAtom._getZ(x)))]
                struct1_atoms.append(atomInd)
                
            struct2_atoms = []
            for x in ringAtoms2:
                x = get_pyatom_from_cppatom(x)
                #Convert from receptor and ligand specific numbering to the combined receptor + ligand numbering
                atomInd = self.close_pdb.origAtomToCombined[str(Point(_StructureAtom._getX(x), _StructureAtom._getY(x),
                                                                      _StructureAtom._getZ(x)))]
                struct2_atoms.append(atomInd)
            
            struct1_res = set([self.close_pdb.atomIDToRes[x] for x in struct1_atoms])
            struct2_res = set([self.close_pdb.atomIDToRes[x] for x in struct2_atoms])
                        
            for res1 in struct1_res:
                if res1 in self.close_pdb.residues:
                    if face_to_face:
                        self.close_pdb.residues[res1].add_pi_stack(struct1_atoms, struct2_atoms)
                    else:
                        self.close_pdb.residues[res1].add_pi_t(struct1_atoms, struct2_atoms)
                    
            for res2 in struct2_res:
                if res2 in self.close_pdb.residues:
                    if face_to_face:
                        self.close_pdb.residues[res2].add_pi_stack(struct1_atoms, struct2_atoms)
                    else:
                        self.close_pdb.residues[res2].add_pi_t(struct1_atoms, struct2_atoms)

        
        '''
        for aromatic1 in self.ligand.aromatics:
            for residue in self.close_pdb.residues.values():
                for aromatic2 in residue.aromatics:
                    if aromatic1.center.dist_to(aromatic2.center) < self.parameters['pi_pi_interacting_dist_cutoff']:
                        residue.add_pi(aromatic2, aromatic1)
        '''

    '''
    Adds pi-cation interactions to the close_pdb() object
    
    Filtered by L2 distance between center of ligand and aromatic ring; projected distance of ligand on ring
    '''
    def get_pi_cation(self):
        '''
        for residue in self.close_pdb.residues.values():
            #Evaluate potential aromatic (from receptor) - cation (from ligand) interactions
            for aromatic in residue.aromatics:
                for atom in self.ligand.atoms:
                    if atom.formal_charge <= 0.1: continue #Cation must be negatively charged
                    if atom.coordinates.dist_to(aromatic.center) < self.parameters['cation_pi_dist_cutoff']:
                        charge_projected = func.project_point_onto_plane(atom.coordinates,aromatic.plane_coeff)
                        if charge_projected.dist_to(aromatic.center) < aromatic.radius + self.parameters['pi_padding_dist']:
                            residue.add_pi_cation(aromatic, atom)

            #Evaluate potential cation (from receptor) - aromatic (from ligand) interactions
            for atom in residue.atoms:
                if atom.formal_charge <= 0.1: continue #Cation must be negatively charged
                for aromatic in self.ligand.aromatics:
                    if atom.coordinates.dist_to(aromatic.center) < self.parameters['cation_pi_dist_cutoff']:
                        charge_projected = func.project_point_onto_plane(atom.coordinates,aromatic.plane_coeff)
                        if charge_projected.dist_to(aromatic.center) < aromatic.radius + self.parameters['pi_padding_dist']:
                            residue.add_pi_cation(aromatic, atom)
        '''        
        picats = interactions.find_pi_cation_interactions(self.close_pdb.receptorOnlySt, struct2=self.ligand.st)
        
        for catAtom, ringAtoms in picats:
            #NOTE: STRUCTURES RETURNED ARE STRUCTURE ATOMS: 
            #http://content.schrodinger.com/Docs/r2017-1/python_api/api/schrodinger.structure._StructureAtom-class.html
            
            struct1_atoms = []
            for x in [catAtom]:
                x = get_pyatom_from_cppatom(x);
                atomInd = self.close_pdb.origAtomToCombined[str(Point(_StructureAtom._getX(x), _StructureAtom._getY(x),
                                                                      _StructureAtom._getZ(x)))]
                struct1_atoms.append(atomInd)
                
            struct2_atoms = []
            for x in ringAtoms:
                x = get_pyatom_from_cppatom(x);
                atomInd = self.close_pdb.origAtomToCombined[str(Point(_StructureAtom._getX(x), _StructureAtom._getY(x),
                                                                      _StructureAtom._getZ(x)))]
                struct2_atoms.append(atomInd)
            
            struct1_res = set([self.close_pdb.atomIDToRes[x] for x in struct1_atoms])
            struct2_res = set([self.close_pdb.atomIDToRes[x] for x in struct2_atoms])
                        
            for res1 in struct1_res:
                if res1 in self.close_pdb.residues:
                    self.close_pdb.residues[res1].add_pi_cation(struct1_atoms, struct2_atoms)
                    
            for res2 in struct2_res:
                if res2 in self.close_pdb.residues:
                    self.close_pdb.residues[res2].add_pi_cation(struct1_atoms, struct2_atoms)
                    
    '''#Included in the contact term! Also, gives crap results.                        
    def get_salt_bridge(self):
        #Schrodinger API lets us use combined atom number, no conversion needed
        for atom1, atom2 in SaltBridgeFinder.find(self.close_pdb.st, group1=self.close_pdb.receptorAtoms,
                                                  group2=self.close_pdb.ligandAtoms):            
            resNum1 = self.close_pdb.atomIDToRes[atom1]
            resNum2 = self.close_pdb.atomIDToRes[atom2]
            
            if resNum1 in self.close_pdb.residues:
                self.close_pdb.residues[resNum1].add_salt_bridge(atom1, atom2)
            
            if resNum2 in self.close_pdb.residues:
                self.close_pdb.residues[resNum2].add_salt_bridge(atom1, atom2)
    '''
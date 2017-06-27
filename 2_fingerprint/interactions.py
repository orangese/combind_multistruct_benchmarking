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

# values taken from Autodock 3.0.5 user guide
# http://autodock.scripps.edu/faqs-help/manual/autodock-3-user-s-guide/AutoDock3.0.5_UserGuide.pdf
# (C_{12} [kcal mol^-1 A^12], C_6 [kcal mol^-1 A^6])
LJ_consts = {
                'C':(2516582.400, 1228.800000),
                'N':(540675.281, 588.245000),
                'O':(230584.301, 429.496730),
                'S':(3355443.200, 1638.400000),
                'H':(81.920, 2.560000)
            }

#e_r_water = 80.1 # at 20 C, wikipedia relative permittivity page
#e_0 = 8.85418782**(-12) # m^-3 kg^-1 s^4 A^2 amps not angstroms!
#electrostatic_const = 1/(4*math.pi*e_r_water*e_0)

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
        self.get_potentials()
        
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
        e_r_water = 80.1 # at 20 C, wikipedia relative permittivity page
        e_0 = 8.85418782*10**(-12) # m^-3 kg^-1 s^4 A^2 amps not angstroms!
        electrostatic_const = 1/(4*math.pi*e_r_water*e_0)
        elem_charge = 1.602*10**(-19)
        angs = 10**(-10)
        kcal_per_mol = (6.022*10**(23)/4184.0)*elem_charge**2 * electrostatic_const/angs
        const = 4.1446 # final answer is in units of kcal/mol, just like LJ
    '''
    def get_potentials(self):
	hydrogens = []
        for lig_atom in self.ligand.all_atoms():
            for residue in self.close_pdb.residues.values():
                for res_atom in residue.atoms:
                    r = lig_atom.dist_to(res_atom)

                    # 1: LJ potential
                    (C12L, C6L) = LJ_consts.get(lig_atom.element,(0,0))
                    (C12R, C6R) = LJ_consts.get(res_atom.element,(0,0))
                    residue.add_lj_potential( (C12L*C12R)**0.5/r**12 - (C6L*C6R)**0.5/r**6 )

                    # 2: electrostatic potential 
                    const = 4.1446
                    q1 = res_atom.formal_charge + res_atom.charge
                    q2 = lig_atom.formal_charge + lig_atom.charge
                    residue.add_electrostatic_potential(q1*q2*const/r)

		    # 3: hydrogen bonds
                    if lig_atom.element in ('N', 'O') and res_atom.element in ('N', 'O'):

                        # find the hydrogen
                        receptor_hydrogens = [atom for atom in list(res_atom.connected_atoms) if atom.element == 'H']
                        ligand_hydrogens = [atom for atom in list(lig_atom.connected_atoms) if atom.element == 'H']
			
                        (best_h, min_dist, donor) = (None, 10.0, None)
                        for h in receptor_hydrogens + ligand_hydrogens:
                            dist = max(lig_atom.dist_to(h), res_atom.dist_to(h))
                            if dist < min_dist:
                                (best_h, min_dist, donor) = (h, dist, h in receptor_hydrogens)
			
			if best_h == None:
                            continue
                        else:
                            hydrogens.append(best_h)
                         
                        angle = math.fabs(180 - func.angle_between_three_points(lig_atom.coordinates,
                                                                                best_h.coordinates, res_atom.coordinates) * 180 / math.pi)
                    
                        score = ( 1/(1+math.exp(4*(min_dist-2.6))) )*( 1/(1+math.exp((angle-60)/10)) )
                        if best_h.score < score:
                           best_h.score = score 
                           best_h.donor = donor
                           best_h.residue = residue
    
        for h in hydrogens:
            h.residue.add_h_bond(h.score,h.donor) 
                        
                        #residue.debug_h(lig_atom,res_atom,best_h,donor,score)


# BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob
# Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite 
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.

# The below is implemented by Joe Paggi and largely based on BINANA

import math_functions as func
import math
from aromatic import Aromatic

class PDB:
    """
    A PDB is resposible for performing actions common to all PDB subclasses.
    These involve assigning bonds and saving entries as PDB files
    """
    protein_resnames = ["ALA", "ARG", "ASN", "ASP", "ASH", "ASX", "CYS", "CYM", "CYX", "GLN", "GLU", "GLH", "GLX", "GLY",
                        "HIS", "HID", "HIE", "HIP", "ILE", "LEU", "LYS", "LYN", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL"]
    def save(self, filename):
        f = open(filename, 'w')
        f.write('\n'.join(atom.CreatePDBLine() for atom in self.all_atoms()))
        f.write('\n')
        f.close()

    def show_hydrophobics(self, filename):
        f = open(filename, 'w')
        f.write('\n'.join(atom.CreatePDBLine() for atom in self.all_atoms() if abs(atom.charge) < .2))
        f.write('\n')
        f.close()

    def show_aromatics(self, filename):
        f = open(filename, 'w')
        f.write('\n'.join(atom.CreatePDBLine() for atom in self.debug_aromatics()))
        f.write('\n')
        f.close()

    def show_interactions(self, filename):
        f = open(filename, 'w')
        f.write('\n'.join(atom.CreatePDBLine() for atom in self.interacting_atoms()))
        f.write('\n')
        f.close()

    # Below are all private
    def _assign_bonds(self):
        for i, atom1 in enumerate(self.all_atoms()[:-1]):
            for atom2 in self.all_atoms()[i+1:]:
                if (atom1.dist_to(atom2) < self._bond_length(atom1.element, atom2.element) * 1.2):
                    atom1.add_neighbor(atom2)
                    atom2.add_neighbor(atom1)

    def _bond_length(self, element1, element2):
        """
        Bond lengths taken from Handbook of Chemistry and Physics. The information provided there was very specific,
        so I tried to pick representative examples and used the bond lengths from those. Sitautions could arise where these
        lengths would be incorrect, probably slight errors (<0.06) in the hundreds.
        """
        element1, element2 = sorted([element1.upper(), element2.upper()])

        UNBONDABLE = ['H', 'CL', 'BR', 'F', 'O', 'P']
        if element1 in UNBONDABLE and element2 in UNBONDABLE: return 0.0

        match = lambda e1, e2: (e1 == element1 and e2 == element2)
        if match('C', 'C'): return 1.53
        if match('N', 'N'): return 1.425
        if match('O', 'O'): return 1.469
        if match('S', 'S'): return 2.048
        if match('C', 'H'): return 1.059
        if match('C', 'N'): return 1.469
        if match('C', 'O'): return 1.413
        if match('C', 'S'): return 1.819
        if match('H', 'N'): return 1.009
        if match('N', 'O'): return 1.463
        if match('O', 'S'): return 1.577
        if match('H', 'O'): return 0.967
        if match('H', 'S'): return 2.025/1.5
        if match('N', 'S'): return 1.633
    
        if match('C', 'F'): return 1.399
        if match('C', 'CL'): return 1.790
        if match('BR', 'C'): return 1.910
        if match('C', 'I'): return 2.162
    
        if match('BR', 'S'): return 2.321
        if match('CL', 'S'): return 2.283
        if match('F', 'S'): return 1.640
        if match('I', 'S'): return 2.687
    
        if match('BR', 'P'): return 2.366
        if match('CL', 'P'): return 2.008
        if match('F', 'P'): return 1.495
        if match('I', 'P'): return 2.490
        if match('O', 'P'): return 1.6
    
        if match('BR', 'N'): return 1.843
        if match('CL', 'N'): return 1.743
        if match('F', 'N'): return 1.406
        if match('I', 'N'): return 2.2
    
        if match('BR', 'SI'): return 2.284
        if match('CL', 'SI'): return 2.072
        if match('F', 'SI'): return 1.636
        if match('P', 'SI'): return 2.264
        if match('S', 'SI'): return 2.145
        if match('SI', 'SI'): return 2.359
        if match('C', 'SI'): return 1.888
        if match('N', 'SI'): return 1.743
        if match('O', 'SI'): return 1.631

        return 1 #If we have no clue, just assume that it's a neighbor, shouldn't affect things?
        assert False, "No bond length info for {} to {}. P.S. I am so sorry, I tried my best.".format(element1, element2)

    def _assign_aromatic_rings(self):
        self.aromatics = []
        PLANAR_ANGLE = 15
        for ring in self._get_rings():
            if any(atom.num_neighbors() == 4 for atom in ring): continue
            is_flat = True
            for t in range(-3, len(ring)-3):
                pt1 = ring[t].coordinates
                pt2 = ring[t+1].coordinates
                pt3 = ring[t+2].coordinates
                pt4 = ring[t+3].coordinates
                # now check the dihedral between the ring atoms to see if it's flat
                angle = func.dihedral(pt1, pt2, pt3, pt4) * 180 / math.pi
                if (angle > -180+PLANAR_ANGLE and angle < -PLANAR_ANGLE) or (angle > PLANAR_ANGLE and angle < 180-PLANAR_ANGLE):
                    is_flat = False
                    break

                # now check the dihedral between the ring atoms and an atom connected to the current atom to see if that's flat too.
                for atom in list(ring[t].connected_atoms):
                    angle = func.dihedral(pt2, pt3, pt4, atom.coordinates) * 180 / math.pi
                    if (angle > PLANAR_ANGLE-180 and angle < -PLANAR_ANGLE) or (PLANAR_ANGLE > 15 and angle < 180-PLANAR_ANGLE):
                        is_flat = False
                        break

            if is_flat: self.aromatics += [Aromatic(list(ring))]

    def _get_rings_containing(self, root):
        """
        Returns all 5 or 6 membered rings containing root
        """
        rings = []
        for i, neighbor1 in enumerate(list(root.connected_atoms)[:-1]):
            for neighbor2 in list(root.connected_atoms)[i+1:]:
                for atom1 in neighbor1.connected_atoms:
                    if atom1 is root: continue
                    for atom2 in neighbor2.connected_atoms:
                        if atom2 is root: continue
                        if atom1 in atom2.connected_atoms: rings += [(root, neighbor1, atom1, atom2, neighbor2)]
                        for atom in list(atom1.connected_atoms):
                            if atom in atom2.connected_atoms:
                                rings += [(root, neighbor1, atom1, atom, atom2, neighbor2)]
        return rings

    def _get_rings(self):
        rings = []
        for root in self.atoms:
            root_rings = self._get_rings_containing(root)
            for ring in root_rings:
                if not any(all(atom in ring2 for atom in ring) for ring2 in rings): rings += [ring]
        return rings

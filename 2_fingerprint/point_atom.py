# BINANA is released under the GNU General Public License (see http://www.gnu.org/licenses/gpl.html).
# If you have any questions, comments, or suggestions, please don't hesitate to contact me, Jacob
# Durrant, at jdurrant [at] ucsd [dot] edu. If you use BINANA in your work, please cite 
# Durrant, J. D. and J. A. McCammon (2011). "BINANA: A novel algorithm for ligand-binding
# characterization." J Mol Graph Model 29(6): 888-893.

# The below is implemented by Joe Paggi and largely based on BINANA

import math
import os
import sys
import textwrap

class Point:
    x=99999.0
    y=99999.0
    z=99999.0
    
    def __init__ (self, x, y ,z):
        self.x = x
        self.y = y
        self.z = z

    def copy_of(self):
        return Point(self.x, self.y, self.z)

    def __str__(self):
        return ("%.3f" % self.x).rjust(8)+("%.3f" % self.y).rjust(8)+("%.3f" % self.z).rjust(8)
        
    def snap(self,reso): # snap the point to a grid
        self.x = round(self.x / reso) * reso
        self.y = round(self.y / reso) * reso
        self.z = round(self.z / reso) * reso
        
    def dist_to(self,apoint):
        return math.sqrt(math.pow(self.x - apoint.x,2) + math.pow(self.y - apoint.y,2) + math.pow(self.z - apoint.z,2))

    def description(self):
        return str(self.x) + " " + str(self.y) + " " + str(self.z)

    def magnitude(self):
        return self.dist_to(Point(0,0,0))

    def CreatePDBLine(self):
        return "HETATM   45      CAU B 400    {}  1.00  0.00      0.0  H".format(str(self))
    def blank():
        #HETATM   45      CAU B 400      28.048   4.763   2.026  1.00  0.00      0.0  H
#        ATOM   TEST    C XXX  24.340   5.139   3.979                       C
        output = "HETATM "
        output = output + str(100).rjust(5) + "C".rjust(5) + "XXX".rjust(4)
        output += '         '
        output = output + ("%.3f" % self.x).rjust(8)
        output = output + ("%.3f" % self.y).rjust(8)
        output = output + ("%.3f" % self.z).rjust(8)
        output = output + "C".rjust(24) 
        return output

class Atom: 
    def __init__ (self):
        self.line = ''
        self.entry = ''
        self.atom_id = -1
        self.atom_name = ''
        self.residue_name = ''
        self.chain = ''
        self.residue_id = -1
        self.coordinates = Point(9999.9, 9999.9, 9999.9)
        self.occupancy = 0.0
        self.b_factor = 0.0
        self.element = ''
        self.formal_charge = 0
        self.charge = 0.0
        self.connected_atoms = set()

    def copy_of(self):
        new = Atom()
        new.line = self.line
        new.entry = self.entry
        new.atom_id = self.atom_id
        new.atom_name = self.atom_name
        new.residue_name = self.residue_name
        new.chain = self.chain
        new.residue_id = self.residue_id
        new.coordinates = self.coordinates.copy_of()
        new.occupancy = self.occupancy
        new.b_factor = self.b_factor
        new.element = self.element
        new.formal_charge = self.formal_charge
        new.charge = self.charge
        return new

    def dist_to(self, other):
        return self.coordinates.dist_to(other.coordinates)

    def string_id(self):
        toreturn = ""
        if self.chain.strip() != '': toreturn = toreturn + self.chain.strip() + ":"
        toreturn = toreturn + self.residue.strip() + '(' + str(self.resid) + "):" + self.atomname.strip() + '(' + str(self.PDBIndex) + ')'
        return toreturn
        
    def CreatePDBLine(self):
        """
        PDB files have a specified number of characters for each field.
        See ftp://ftp.wwpdb.org/pub/pdb/doc/format_descriptions/Format_v33_Letter.pdf
        for a complete description.
        Left vs right justified is consistent with Maestro produced files

        Charge occupies an otherwise blank part of the line. Is not part of PDB format.
        """
        output = self.entry.ljust(6)
        output += str(self.atom_id).rjust(5)
        output += ' '
        if len(self.atom_name) == 4:
            output += self.atom_name
        else:
            output += ' ' + self.atom_name.ljust(3)
        output += ' '
        output += self.residue_name.rjust(3)
        output += ' '
        output += self.chain.rjust(1)
        output += str(self.residue_id).rjust(4)
        output += ' ' * 4
        output += str(self.coordinates)
        output += self.occupancy
        output += self.b_factor
        output += ' '
        output += str(self.charge)[:8].rjust(8)
        output += ' '
        output += self.element.rjust(2)
        return output

    def num_neighbors(self):
        return len(self.connected_atoms)

    def add_neighbor(self, atom):
        self.connected_atoms.add(atom)
    
    def SideChainOrBackBone(self):
        return 'BACKBONE' if self.atomname in ('CA', 'C', 'O', 'N') else 'SIDECHAIN'

    def from_schrod(self, atom):
        self.line = ''
        self.entry = 'ATOM'
        self.atom_id = atom.index
        self.atom_name = atom.pdbname.strip()
        self.residue_name = atom.pdbres.strip()
        self.chain = atom.chain.strip()
        self.residue_id = atom.resnum
        self.coordinates = Point(atom.x, atom.y, atom.z)
        self.occupancy = ' ' * 6
        self.b_factor  = ' ' * 6
        self.element = atom.element.strip()
        self.formal_charge = atom.formal_charge
        self.charge = atom.partial_charge
    
    def ReadPDBLine(self, line):
        self.line = line
        self.entry = line[:6].strip()
        self.atom_id = int(line[6:11].strip())
        self.atom_name = line[12:16].strip()
        self.residue_name = line[17:20].strip()[-3:]
        self.chain = line[21]
        self.residue_id = int(line[22:26].strip())
        self.coordinates = Point(float(line[30:38]), float(line[38:46]), float(line[46:54]))
        self.occupancy = line[54:60]
        self.b_factor = line[60:66]
        self.element = line[76:78].strip()
        formal_charge = line.strip()[78:]
        if formal_charge:
            if formal_charge[1] == '-':
                self.formal_charge = - int(formal_charge[0])
            else:
                self.formal_charge = int(formal_charge[0])
        else:
            self.formal_charge = 0
        self.charge = float(self.formal_charge)

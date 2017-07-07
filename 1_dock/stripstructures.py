#!/share/PI/rondror/software/schrodinger2017-1/run
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
from schrodinger.structure import _StructureAtom
from schrodinger.structure import Structure
from schrodinger.structutils.structalign import StructAlign
import os
from os import listdir
from os.path import isfile, join
import sys

def strip():
    dirFiles = [f for f in listdir(".") if isfile(join("./", f)) and (f.endswith(".pdb") or f.endswith(".mae"))]
    dirPath = os.getcwd()
    structs = []
    
    print("--Reading in structs...")
    for f in dirFiles:
        struct = StructureReader("./" + f).next()

        #Delete all atoms that do not belong to chain A, removes multiple copies of the protein-ligand complex
        nonAAtoms = [] 
        for atom in struct.atom:
            if atom._getAtomChain() != "A":
                nonAAtoms.append(atom.__int__())
        struct.deleteAtoms(nonAAtoms)
        structs.append(struct)

    print("--Removing waters from structs...")
    newStructs = [] #List for structs without water
    for struct in structs:
        delAtomList = []
        for mol in struct.molecule:
            atomList = map(lambda x: x.element, mol.atom)
            if(len(atomList) == 1 and atomList[0] == 'O'): #Found a water molecule
                delAtomList.append(_StructureAtom.__int__(mol.atom[0]))
        struct.deleteAtoms(delAtomList)
        newStructs.append(struct)

    print("--Aligning structures...")
    structAlign = StructAlign()
    structAlign.align(newStructs[0], newStructs[1:len(structs)])
    
    for struct in newStructs:
        st_writer = StructureWriter("./" + struct._getTitle() + ".mae")
        st_writer.append(struct)
        st_writer.close()
 

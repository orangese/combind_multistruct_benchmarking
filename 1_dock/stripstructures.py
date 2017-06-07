#!/share/PI/rondror/software/schrodinger2016-1/run
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
    os.chdir("./raw_pdbs")
    dirFiles = [f for f in listdir(".") if isfile(join("./", f)) and (f.endswith(".pdb") or f.endswith(".mae"))]
    dirPath = os.getcwd()
    structs = []
    
    print("--Reading in structs...")
    for f in dirFiles:
        struct = StructureReader("./" + f).next()
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
    
    for ind, struct in enumerate(newStructs):
        if struct._getTitle() is not "":
            st_writer = StructureWriter("./" + struct._getTitle() + ".mae")
        else:
            st_writer = StructureWriter("./" + os.path.splitext(dirFiles[ind])[0].upper() + ".mae")
        st_writer.append(struct)
        st_writer.close()



   # print("--Aligning structs...")
   # structAlign = StructAlign()
   # structAlign.align(newStructs[0], newStructs[1:len(structs)])
    
   # for ind, struct in enumerate(newStructs):
   #     if struct._getTitle() is not "":
   #         st_writer = StructureWriter("./" + struct._getTitle() + ".mae")
   #    else:
   #         st_writer = StructureWriter("./" + os.path.splitext(dirFiles[ind])[0].upper() + ".mae")
   #     st_writer.append(struct)
   #     st_writer.close()
    
    os.system("mkdir ../stripped")
    os.system("mv *.mae ../stripped")
    os.chdir("../")

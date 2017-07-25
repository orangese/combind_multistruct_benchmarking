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
    rawFiles = [f for f in listdir("raw_pdbs") if f.endswith(".pdb") and not os.path.exists('stripped/{}.mae'.format(f.split('.')[0]))]
    print rawFiles
    structs = []
     
    print("--Reading in structs...")
    for f in rawFiles:
        struct = StructureReader("raw_pdbs/" + f).next()
        
        #If there are multiple chains in our PDB File, delete all non-A chains
        chainNames = set([chain._getChainName() for chain in struct.chain])
        if len(chainNames) > 1 and "A" in chainNames:
            nonAAtoms = [] 
            for atom in struct.atom:
                if atom._getAtomChain() not in ['A', 'L', 'Z']:# the ligand is sometimes in L or Z
                    nonAAtoms.append(atom.__int__())
            struct.deleteAtoms(nonAAtoms)

        #If the PDB lacks a title, give it the same title as the filename    
        if struct._getTitle() == "" or struct._getTitle() == 'xxxx':
            struct._setTitle(os.path.splitext(f)[0].upper())
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

    if 'stripped' in listdir('.') and len(listdir('stripped')) > 0:
        template = StructureReader('stripped/{}'.format(listdir('stripped')[0])).next()
        structAlign.align(template, newStructs)
    else:
        structAlign.align(newStructs[0], newStructs[1:])
    
    os.system('mkdir -p stripped')
    for struct in newStructs:
        st_writer = StructureWriter("./stripped/" + struct._getTitle() + ".mae")
        st_writer.append(struct)
        st_writer.close()


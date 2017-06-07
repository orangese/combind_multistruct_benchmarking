import sys
import os
import multiprocessing as mp
import xlrd
from sets import Set
import wget
import ssl
import multiprocessing as mp
from tqdm import tqdm

SCHRODINGER = "/share/PI/rondror/software/schrodinger2016-1"

def sdToMae(f):
    os.system(SCHRODINGER + "/ligprep -WAIT -adjust_itc -epik -isd " + f + " -omae " + f[:-3] + ".mae")

def parseExcelMapFile(inputFile, parseSmiles):
    numCores = 4
    ssl._create_default_https_context = ssl._create_unverified_context

    os.chdir("./raw_pdbs")
    sh = xlrd.open_workbook(inputFile).sheet_by_index(0)

    pdbList = Set([])
    smilesList = []
    smilesIDList = []

    pdbColumnIndex = -1
    smilesColumnIndex = -1
    smileIDIndex = -1
    
    #Find the correct column numbers
    for colnum in range(sh.ncols):
        if(sh.cell(0, colnum).value == "PDB ID"):
            pdbColumnIndex = colnum
        if(sh.cell(0, colnum).value == "SMILES"):
            smilesColumnIndex = colnum
        if("Compound" in sh.cell(0, colnum).value and "ID" in sh.cell(0,colnum).value):
            smilesIDIndex = colnum

    #Iterate through and collect PDB IDs
    for rownum in range(sh.nrows):
        tempValue = sh.cell(rownum, pdbColumnIndex).value
        if tempValue != "PDB ID" and tempValue != "":
            pdbList.add(tempValue)

    for rownum in range(sh.nrows):
        tempValue = sh.cell(rownum, smilesColumnIndex).value
        if tempValue != "SMILES" and tempValue != "":
            smilesList.append(tempValue)
    
    for rownum in range(sh.nrows):
        tempValue = sh.cell(rownum, smilesIDIndex).value
        if tempValue != "Compound_ID" and tempValue != "Compound ID" and tempValue != "":
            smilesIDList.append(tempValue)
    
    if(parseSmiles):
        os.chdir("..")
        os.mkdir("exp_ligands")
        os.chdir("exp_ligands")

        pool = mp.Pool(numCores)
        print("Converting SMILE Strings...")

        for _ in tqdm(pool.imap_unordered(convertSmileFile, zip(smilesList, smilesIDList)), total=len(zip(smilesList, smilesIDList))):
            pass

        os.system("rm *.smi *.log *.tar.gz")

        os.chdir("../raw_pdbs")
    
    with open("pdb_list.txt", "w") as f:
        for pdbID in pdbList:
            f.write(pdbID + '\n')
    f.close()
    os.chdir("../")
    #parseTextMapFile("pdb_list.txt")

def convertSmileFile(smileData):
    smileString = smileData[0]
    smileID = smileData[1]
    
    with open(smileID + ".smi", "w") as f:
        f.write(smileString + " \'" + smileID+"\'")
        newSmileString = smileString.replace("_","") #Remove all instances of "_" since the script looks for "_ligand" 
    os.system("/share/PI/rondror/software/schrodinger2016-1/ligprep -WAIT -adjust_itc -epik -ismi " + smileID + ".smi -omae " + newSmileString + "_ligand.mae")

def parseTextMapFile(inputFile):
    os.chdir("raw_pdbs")
    with open(inputFile, 'r') as f:
        for line in f:
            line = str.strip(line)
            tempDownload = wget.download("http://files.rcsb.org/download/" + line + ".pdb", out=".")
    os.chdir("../")

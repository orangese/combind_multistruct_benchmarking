#!/share/PI/rondror/software/schrodinger2017-1/run
print("!Importing necessary Schrodinger Libraries!")
import sys
import os

import parsemap
import stripstructures
import ligprep #method extract is thread safe
import gridgen #method runGlide is thread safe
import processing #method process is thread safe
import dock

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"
DATA = "/scratch/PI/rondror/docking_data"

os.chdir(DATA)

commands = sys.argv[1]
structures = sys.argv[2:]

toRun = []
if("a" in commands):
    toRun = ["r","s","p","l","g","x"] #Raw PDBs, Strip Structures, Process Structures, Grids, Dock
else:
    toRun = list(commands[1:])

for struct in structures:#Go through the given structures, performing the commands specified
    
    if(not os.path.exists("./" + struct)): #Make sure that the inputted directory exists
        raise Exception("Error: ./" + struct + " does not exist!")

    os.chdir("./"+struct)#Note: All work will be done with the current working directory in the current struct folder

    #Download "raw" structures from the PDB
    if("r" in toRun):
        print("!Downloading Raw PDB Files")
        #Find potential files that could contain the structures
        structFiles = [name for name in os.listdir(".") if os.path.isfile(name)]
        structFiles = list(filter(lambda x: x.endswith(".txt"), structFiles))

        #Use the first .txt file we find to download structures from
        os.system("mkdir raw_pdbs")
        os.system("mv " + structFiles[0] + " raw_pdbs")
        fileExtension = os.path.splitext(structFiles[0])[1]
        parsemap.parseTextMapFile(structFiles[0])

    #Strip "raw" structures of waters and align them (not parallelized since alignment requires all structures)
    if("s" in toRun):
        print("!Stripping Raw PDB Files")
	stripstructures.strip()
	
    if("p" in toRun):
        print("!Processing Structures")
        processing.process()

    if("l" in toRun):#Extract Ligands, assumes processed folder exists
        print("!Extracting Ligands")
        ligprep.extractLigands()

    if("g" in toRun):#Generate grids, assumes processed folder exists
        print("!Generating Grids")
        if("processed" not in [name for name in os.listdir(".") if os.path.isdir(name)]):
            raise Exception("Error: -g (generate grids) requires the folder processed to be the parent structure directory")

        os.system("mkdir -p grids")
        #makeGrids = [x.split('.')[0] for x in os.listdir('processed') if not os.path.exists(os.getcwd()+'/grids/'+x.split('.')[0]+'/'+x.split('.')[0]+'.zip')]
        makeGrids = ['4OBP','4OBQ','4U43','4U45']
    
        for i in makeGrids:
            os.system('cp ./processed/'+i+'.mae grids')
        
	os.chdir("./grids")

        gridgen.generateIn(makeGrids)
        results = gridgen.runGlides([x+'.in' for x in makeGrids])
        
        os.system("rm *.log *.mae *.in gpu* *.json")

        for x in os.listdir('.'):
            if x[-4:] == '.zip':
                os.system("mkdir -p " + x[:-4])
                os.system("mv " + x + " ./" + x[:-4] + "/")

        os.chdir("..")

    if("x" in toRun):
        os.system('mkdir -p xglide')
        dock.dockDataset(struct)


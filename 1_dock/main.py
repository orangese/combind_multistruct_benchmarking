#!/share/PI/rondror/software/schrodinger2017-1/run
print("!Importing necessary Schrodinger Libraries!")
import sys
import os
import multiprocessing as mp

#Self-Defined Helper Libraries, Note: All functions in these libraries assume that you are in the parent directory for the given struct (ex. /docking/data/A2AR/)
import parsemap
import stripstructures
import ligprep #method extract is thread safe
import gridgen #method runGlide is thread safe
import processed #method process is thread safe

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"
GLIDE = "/share/PI/rondror/software/schrodinger2017-1/glide"

DATA = "/scratch/PI/rondror/docking_data"
DOCKING_SCRIPT = "/share/PI/rondror/$USER/combind/1_dock/dock_ligand_dir_to_grid_dir.sh" 

os.chdir(DATA)

commands = sys.argv[1]
structures = sys.argv[2:]
numCores = 4
print("!Initializing: Found " + str(4) + " cores!")
pool = mp.Pool(4)

toRun = []
if("a" in commands):
    toRun = ["r","s","p","l","g","d"] #Raw PDBs, Strip Structures, Process Structures, Grids, Dock
else:
    toRun = list(commands[1:])

for struct in structures:#Go through the given structures, performing the commands specified
    
    if(not os.path.exists("./" + struct)): #Make sure that the inputted directory exists
        raise Exception("Error: ./" + struct + " does not exist!")

    os.chdir("./"+struct)#Note: All work will be done with the current working directory in the current struct folder

    if("r" in toRun):#Download Raw PDBs
        print("!Downloading Raw PDB Files")
        structFiles = ["raw_pdbs.txt"]
        os.system("mkdir raw_pdbs")
        os.system("mv " + structFiles[0] + " raw_pdbs")
        fileExtension = os.path.splitext(structFiles[0])[1]
        if(fileExtension == ".xls" or fileExtension == ".xlsx"):
            parsemap.parseExcelMapFile(structFiles[0], True)#True = parseSmiles
        else:
            parsemap.parseTextMapFile(structFiles[0])

    if("s" in toRun):#Stripped assumes that there exists a folder called "raw_pdbs" in the parent struct directory
        print("!Stripping Raw PDB Files")
        currentDirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
        if("raw_pdbs" not in currentDirectories):
            raise Exception("Error: -s (strip structures) requires the folder raw_pdbs to be the parent structure directory")
        os.system('mkdir stripped')
	os.system('cp raw_pdbs/*.pdb stripped/')
	os.chdir('stripped')
	stripstructures.strip()
	os.system('rm *.pdb')
	os.chdir('..')

    if("p" in toRun):#Process Structures assumes that there exists a "stripped" folder in the parent struct directory
        print("!Processing Stripped Files")
        currentDirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
        if("stripped" not in currentDirectories):
            raise Exception("Error: -p (process structures) requires the folder stripped to be the parent structure directory")
        os.system("mkdir processed")
	os.system("cp stripped/*.mae processed/")
	os.chdir("./processed")
        pool = mp.Pool(numCores)
        processed.process(pool) # this function (1) assumes you're in the processed dir (2) copies files from stripped
        os.system("rm temp*")
	os.chdir("..")

    if("l" in toRun):#Extract Ligands, assumes processed folder exists
        print("!Extracting Ligands")
        currentDirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
        if("processed" not in currentDirectories):
            raise Exception("Error: -l (extract ligands) requires the folder processed to be the parent structure directory")
        os.system("mkdir ligands")
        processedFiles = [f for f in os.listdir("./processed") if os.path.isfile(os.path.join("./processed", f))]
        processedFiles = map(lambda x: os.path.splitext(x)[0], processedFiles)
        pool = mp.Pool(numCores)
        print(processedFiles)
        results = pool.map(ligprep.extract, processedFiles)

    if("g" in toRun):#Generate grids, assumes processed folder exists
        print("!Generating Grids")
        currentDirectories = [name for name in os.listdir(".") if os.path.isdir(name)]
        if("processed" not in currentDirectories):
            raise Exception("Error: -g (generate grids) requires the folder processed to be the parent structure directory")
        os.system("mkdir grids")
        os.system("cp ./processed/* grids")
        currentInFiles = [f for f in os.listdir("./grids") if os.path.isfile(os.path.join("./grids", f))]
        currentInFiles = map(lambda x: os.path.splitext(x)[0], currentInFiles)
	os.chdir("./grids")
        pool = mp.Pool(numCores)
        gridgen.generateIn(currentInFiles,pool)
        currentInFiles = map(lambda x: x + ".in", currentInFiles)
        results = pool.map(gridgen.runGlide, currentInFiles)

        os.system("rm *.log")
        os.system("rm *.mae")
        os.system("rm *.in")

        files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

        for f in files:
            os.system("mkdir " + str(f[:-4]))
            os.system("mv " + f + " ./" + f[:-4] + "/" + f)

        os.chdir("..")

    if("d" in toRun):#Submit the docking run
        print("!Submitting docking run")
        print(os.getcwd())
        os.system('mkdir glide')
        os.system(DOCKING_SCRIPT + " " + os.getcwd()+"/grids " + os.getcwd() +"/ligands " + os.getcwd() + "/glide/")

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
import processing #method process is thread safe

SCHRODINGER = "/share/PI/rondror/software/schrodinger2017-1"
GLIDE = "/share/PI/rondror/software/schrodinger2017-1/glide"

HERE = os.getcwd() + '/'

DATA = "/scratch/PI/rondror/docking_data"
DOCKING_SCRIPT = HERE + "dock_ligand_dir_to_grid_dir.sh" 
PROCESSING_SCRIPT = HERE + "processing_script.sh"
RMSD_SCRIPT= HERE + 'compute_rmsds.py'
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
        currentDirectories = [name for name in os.listdir(".") if os.path.isdir(name)]

        #Stripped assumes that there exists a folder called "raw_pdbs" in the parent struct directory
        if("raw_pdbs" not in currentDirectories):
            raise Exception("Error: -s (strip structures) requires the folder raw_pdbs to be the parent structure directory")

        os.system('mkdir stripped')
	    os.system('cp raw_pdbs/*.pdb stripped/')
	    os.chdir('stripped')
	    stripstructures.strip()
	    os.system('rm *.pdb')
	    os.chdir('..')

    if("p" in toRun):
        print("!Processing Structures")
        processing.process()

    if("l" in toRun):#Extract Ligands, assumes processed folder exists
        print("!Extracting Ligands")

        if("processed" not in  [name for name in os.listdir(".") if os.path.isdir(name)]):
            raise Exception("Error: -l (extract ligands) requires the folder processed to be the parent structure directory")

        os.system("mkdir ligands")
        processedFiles = [f for f in os.listdir("./processed") if os.path.isfile(os.path.join("./processed", f))]
        processedFiles = map(lambda x: os.path.splitext(x)[0], processedFiles)
        #Extracting ligands is very lightweight so this is single-threaded
        ligprep.extractLigands(processedFiles)

    if("g" in toRun):#Generate grids, assumes processed folder exists
        print("!Generating Grids")
        if("processed" not in [name for name in os.listdir(".") if os.path.isdir(name)]):
            raise Exception("Error: -g (generate grids) requires the folder processed to be the parent structure directory")

        os.system("mkdir grids")
        os.system("cp ./processed/* grids")
        currentInFiles = [f for f in os.listdir("./grids") if os.path.isfile(os.path.join("./grids", f))]
        currentInFiles = map(lambda x: os.path.splitext(x)[0], currentInFiles)
	    os.chdir("./grids")

        gridgen.generateIn(currentInFiles, pool)
        currentInFiles = map(lambda x: x + ".in", currentInFiles)
        results = gridgen.runGlides(currentInFiles)

        os.system("rm *.log")
        os.system("rm *.mae")
        os.system("rm *.in")

        files = [f for f in os.listdir(".") if os.path.isfile(os.path.join(".", f))]

        for f in files:
            os.system("mkdir " + str(f[:-4]))
            os.system("mv " + f + " ./" + f[:-4] + "/" + f)

        os.chdir("..")

    if("d" in toRun):#Submit the docking run
        missing_glide = 0
        total_glide = 0
        os.system('mkdir -p glide')

        glidesExist = map(lambda x: os.path.exists(os.getcwd()+'/glide/'+x+'/'+x+'_pv.maegz'), os.listdir('glide'))

        if len(os.listdir('glide')) > 0 and glidesExist.count(False) == 0:
            print 'Docking results already exist.'
        else:
            if missing_glide > 0: 
                print 'Missing ' + str(missing_glide) + ' of ' + str(total_glide) + ' docking results.'
            print 'Submitting docking jobs...'
            os.system(DOCKING_SCRIPT + " " + os.getcwd()+"/grids " + os.getcwd() +"/ligands " + os.getcwd() + "/glide/")

    if('m' in toRun):
        os.chdir(HERE)
        os.system('sbatch --time=2:00:00 --job-name=r-'+struct+' -n 1 -p rondror '+RMSD_SCRIPT+' '+struct+' '+HERE)


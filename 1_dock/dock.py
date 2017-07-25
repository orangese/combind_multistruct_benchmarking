#!/share/PI/rondror/software/schrodinger2016-1/run
import os
import slurm
from multiprocessing import Pool

SCHRODINGER = '/share/PI/rondror/software/schrodinger2017-1'
DATA_DIR = '/scratch/PI/rondror/docking_data'
LIGANDS_DIR = 'ligands'
PROCESSED_DIR = 'processed'
GLIDE_DIR = 'xglide'
GRIDS_DIR = 'grids'

XGLIDE_IN = '''GRIDFILE   {}_grid.zip
LIGANDFILE   {}_ligand.mae
USE_REF_LIGAND   True
REF_LIGAND_FILE   {}_ligand.mae
CORE_DEFINITION   allheavy
DOCKING_METHOD   confgen
EXPANDED_SAMPLING   True
NENHANCED_SAMPLING   4
POSES_PER_LIG   300
POSTDOCK_NPOSE   300
NMAXRMSSYM   300
MAXKEEP   100000
MAXREF   1000
PRECISION   SP
WRITEREPT   True
'''

def glideExists(dataset, ligand, grid):
    return os.path.exists("{}/{}/{}/{}_ligand-to-{}/{}_ligand-to-{}_pv.maegz".format(
        DATA_DIR, dataset, GLIDE_DIR, ligand, grid, ligand, grid)) or os.path.exists("{}/{}/{}/{}_ligand-to-{}/{}_ligand-to-{}-out.maegz".format(
                    DATA_DIR, dataset, GLIDE_DIR, ligand, grid, ligand, grid))

def glideFailed(dataset, ligand, grid):
    ligToGrid = "{}_ligand-to-{}".format(ligand, grid)
    ligGridDockDir = "{}/{}/{}/{}".format(DATA_DIR, dataset, GLIDE_DIR, ligToGrid)
    logFile = "{}/{}.log".format(ligGridDockDir, ligToGrid)

    if os.path.exists(logFile):
        if "Total elapsed time" in open(logFile).read(): #If we got to the end of Glide, something messed up
            return True

def dock(dataset, ligand, grid):
    #Setup the necessary files for docking, remove old ones if needed
    ligToGrid = "{}_ligand-to-{}".format(ligand, grid)
    ligGridDockDir = "{}/{}/{}/{}".format(DATA_DIR, dataset, GLIDE_DIR, ligToGrid)

    if os.path.exists(ligGridDockDir):
        os.system("rm -rf {}".format(ligGridDockDir))

    os.system("mkdir {}".format(ligGridDockDir))

    gridFile = "{}/{}/{}/{}/{}.zip".format(DATA_DIR, dataset,  GRIDS_DIR, grid, grid)
    os.system("cp {} {}/{}_grid.zip".format(gridFile, ligGridDockDir, ligToGrid))

    ligandFile = "{}/{}/{}/{}_ligand.mae".format(DATA_DIR, dataset, LIGANDS_DIR, ligand)
    os.system("cp {} {}/{}_ligand.mae".format(ligandFile, ligGridDockDir, ligToGrid))
    
    #Write out the input file
    with open("{}/{}.in".format(ligGridDockDir, ligToGrid), 'w+') as f:
        f.write(XGLIDE_IN.format(ligToGrid, ligToGrid, ligToGrid))

    #Submit the docking job and wait for it to finish
    os.chdir(ligGridDockDir)
    slurm.salloc("{}/glide {}/{}.in -WAIT".format(SCHRODINGER, ligGridDockDir, ligToGrid), 
            1, "02:30:00") #jobString, numProcessors, timeLimit

    #Return (ligand,grid) so we know the directory to check for success
    return (ligand, grid)

def dockHelper(ligGrid):
    return dock(*ligGrid)

def dockDataset(dataset):
    gridsDir = "{}/{}/{}".format(DATA_DIR, dataset, GRIDS_DIR)
    ligandsDir = "{}/{}/{}".format(DATA_DIR, dataset, LIGANDS_DIR)

    #Contain the structure name (without _ligand or extensions)
    structures = [o for o in os.listdir(gridsDir) if os.path.isdir(os.path.join(gridsDir,o))]

    ligands = [os.path.splitext(o)[0].replace("_ligand","") for o in os.listdir(ligandsDir) if os.path.isfile(os.path.join(ligandsDir,o))]
    ligands = [l for l in ligands if l not in structures]# do not dock ligands that we have structures for
        
    toDock = [] #List of (ligand, grid) tuples that need to be submitted for docking

    for ligand in ligands:
        for grid in structures:
            if not glideExists(dataset, ligand, grid):
                toDock.append((dataset, ligand, grid))

    print("{} ligand-grid pairs are missing! Checking for failures...".format(str(len(toDock))))
    
    dockFailures = [] #List of structures that failed because docking couldn't find a good pose

    for dataset, ligand, grid in toDock:
        if glideFailed(dataset, ligand, grid):
            dockFailures.append((dataset, ligand, grid))

    print("Of {} missing ligand-grid pairs, {} were failures (listed below). Not submitting these...".format(str(len(toDock)), str(len(dockFailures))))
    print([(x[1], x[2]) for x in dockFailures])

    toDock = [x for x in toDock if x not in dockFailures]
    print("Submitting:")
    print([(x[1], x[2]) for x in toDock])

    pool = Pool(6) #We have 6 Glide SP Docking licenses

    while len(toDock) != 0:
        currentlyDocking = toDock
        toDock = []

        for finishedLigand, finishedGrid in pool.imap_unordered(dockHelper, currentlyDocking):
            #Check to make sure that glide worked successfully
            if not glideExists(dataset, finishedLigand, finishedGrid) and not glideFailed(dataset, finishedLigand, finishedGrid):
                toDock.append((dataset, finishedLigand, finishedGrid))

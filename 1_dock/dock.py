#!/share/PI/rondror/software/schrodinger2016-1/run
import os
import slurm
from multiprocessing import Pool

SCHRODINGER = '/share/PI/rondror/software/schrodinger2017-1'
DATA_DIR = '/scratch/PI/rondror/docking_data'
LIGANDS_DIR = 'ligands'
GRIDS_DIR = 'grids'

XGLIDE_IN = '''GRIDFILE   {}_grid.zip
LIGANDFILE  {}_ligand.mae
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

GLIDE_IN = '''GRIDFILE   {}_grid.zip
LIGANDFILE  {}_ligand.mae
NENHANCED_SAMPLING   4
POSES_PER_LIG   150
POSTDOCK_NPOSE   150
PRECISION   SP
WRITEREPT   True
'''

def glideExists(dataset, ligand, grid):
    return os.path.exists("{}/{}/{}/{}_ligand-to-{}/{}_ligand-to-{}_pv.maegz".format(
        DATA_DIR, dataset, GLIDE_DIR, ligand, grid, ligand, grid))

def dock(dataset, ligand, grid, xDock = True): #xDock = Extra Effort Docking
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

    with open("{}/{}.in".format(ligGridDockDir, ligToGrid), 'w+') as f:
        if xDock:
            f.write(XGLIDE_IN.format(ligToGrid, ligToGrid))
        else:
            f.write(GLIDE_IN.format(ligToGrid, ligToGrid))

    #Submit the docking job and wait for it to finish
    os.chdir(ligGridDockDir)
    slurm.salloc("{}/glide {}/{}.in -WAIT".format(SCHRODINGER, ligGridDockDir, ligToGrid), 
            1, "10:00:00") #jobString, numProcessors, timeLimit

    #Return (ligand,grid) so we know the directory to check for success
    return (ligand, grid)

def dockHelper(ligGrid):
    return dock(*ligGrid)

def dockDataset(dataset, xDock=True):
    print("Docking {}...".format(dataset))
    if xDock:
        GLIDE_DIR = 'xglide'
    else:
        GLIDE_DIR = 'glide'

    gridsDir = "{}/{}/{}".format(DATA_DIR, dataset, GRIDS_DIR)

    structures = [o for o in os.listdir(gridsDir) if os.path.isdir(os.path.join(gridsDir,o))]
    
    toDock = [] #List of (ligand, grid) tuples that need to be submitted for docking

    for ligand in structures:
        for grid in structures:
            if not glideExists(dataset, ligand, grid):
                toDock.append((dataset, ligand, grid, xDock))

    print("{} ligand-grid pairs are missing! Submitting docking...".format(str(len(toDock))))

    pool = Pool(8) #We have about 8 glide licenses

    while len(toDock) !=0:
        currentlyDocking = toDock
        toDock = []

        for finishedLigand, finishedGrid in pool.imap_unordered(dockHelper, currentlyDocking):
            #Check to make sure that glide worked successfully
            if not glideExists(dataset, finishedLigand, finishedGrid):
                toDock.append((finishedLigand, finishedGrid, xDock))
                print("{} to {} failed! Resubmitting to the queue...".format(finishedLigand, finishedGrid))

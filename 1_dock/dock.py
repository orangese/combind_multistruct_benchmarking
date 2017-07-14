#!/share/PI/rondror/software/schrodinger2016-1/run
import os
import slurm
from multiprocessing import Pool
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher

SCHRODINGER = '/share/PI/rondror/software/schrodinger2017-1'
DATA_DIR = '/scratch/PI/rondror/docking_data'
LIGANDS_DIR = 'ligands'
PROCESSED_DIR = 'processed'
GLIDE_DIR = '' #To be determined by xDock
GRIDS_DIR = 'grids'

#.format(receptor, resName, resName, resName, ligToGrid, ligToGrid)
IFGLIDE_IN = '''INPUT_FILE {}.mae

STAGE VDW_SCALING
  BINDING_SITE ligand {}

STAGE PREDICT_FLEXIBILITY
  BINDING_SITE ligand {}

STAGE INITIAL_DOCKING
  BINDING_SITE ligand {}
  INNERBOX 10.0
  OUTERBOX auto
  LIGAND_FILE {}_ligand.mae
  LIGANDS_TO_DOCK all
  VARIANTS_TO_RUN A,B,C,D,E,F,G
  DOCKING_RINGCONFCUT 2.5
  DOCKING_AMIDE_MODE penal

STAGE COMPILE_RESIDUE_LIST
  DISTANCE_CUTOFF 5.0

STAGE PRIME_REFINEMENT
  NUMBER_OF_PASSES  1
  USE_MEMBRANE no
  OPLS_VERSION OPLS_2005

STAGE GLIDE_DOCKING2
  BINDING_SITE ligand Z:999
  INNERBOX 5.0
  OUTERBOX auto
  LIGAND_FILE  {}_ligand.mae
  LIGANDS_TO_DOCK existing
  DOCKING_PRECISION SP
  DOCKING_CANONICALIZE False
  DOCKING_RINGCONFCUT 2.5
  DOCKING_AMIDE_MODE penal

STAGE SCORING
  SCORE_NAME  r_psp_IFDScore
  TERM 1.000,r_psp_Prime_Energy,1
  TERM 9.057,r_i_glide_gscore,0
  TERM 1.428,r_i_glide_ecoul,0
  REPORT_FILE report.csv
'''

XGLIDE_IN = '''GRIDFILE   {}_grid.zip
LIGANDFILE  {}_ligand.mae
DOCKING_METHOD confgen
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

def glideFailed(dataset, ligand, grid):
    ligToGrid = "{}_ligand-to-{}".format(ligand, grid)
    ligGridDockDir = "{}/{}/{}/{}".format(DATA_DIR, dataset, GLIDE_DIR, ligToGrid)
    logFile = "{}/{}.log".format(ligGridDockDir, ligToGrid)

    if os.path.exists(logFile):
        if "Total elapsed time" in open(logFile).read(): #If we got to the end of Glide, something messed up
            return True

def dock(dataset, ligand, grid, xDock = True, inducedFit = False): #xDock = Extra Effort Docking, inducedFit = Induced Fit Docking; Must be either both False or One True; Otherwise default
    print((xDock, inducedFit))
    #Make sure that ONLY ONE docking method has been chosen; Defaults to regular docking
    assert ((not xDock and not inducedFit) or (xDock != inducedFit))

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
   
    #Extra steps are necessary for induced fit docking
    if inducedFit:
        #Induced fit docking requires an extra file to be in the directory
        processedFile = "{}/{}/{}/{}.mae".format(DATA_DIR, dataset, PROCESSED_DIR, grid)
        os.system("cp {} {}/{}.mae".format(processedFile, ligGridDockDir, grid))
        
        #The Induced Fit docking input file requires the chain name and residue number of the ligand
        #NOTE: Here we assume that the ligand has one chain and one residue, this should be the case
        processedStructure = StructureReader(processedFile).next() #Read in processed structure

        # Identify signal ligand, we need the ligand's location IN the receptor, can't use ligand .mae file
        ligand = None

        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(processedStructure)

        assert len(ligands) != 0, 'Error: Could not find a ligand for {}'.format(struct_base)

        if(len(ligands) > 1):
            print('Warning: There are multiple ligand sized molecules, picking the first one')

        ligandStructure = ligands[0].st

        ligandAtom = ligandStructure.atom[1] #Get the first (index starts at 1) atom, find its residue num and chain name
        ligResidueNum = ligandAtom.getResidue()._getResnum()
        ligChainName = ligandAtom.getChain()._getChainName()
        resNum = str(ligChainName)+":"+str(ligResidueNum)
    
    #Write out the input file
    if inducedFit:
        with open("{}/{}.inp".format(ligGridDockDir, ligToGrid), 'w+') as f: #Induced Fit docking requires a .inp file
            f.write(IFGLIDE_IN.format(grid, resNum, resNum, resNum, ligToGrid, ligToGrid))
    else:
        with open("{}/{}.in".format(ligGridDockDir, ligToGrid), 'w+') as f:
            if xDock:
                f.write(XGLIDE_IN.format(ligToGrid, ligToGrid))
            else:
                f.write(GLIDE_IN.format(ligToGrid, ligToGrid))

    #Submit the docking job and wait for it to finish
    os.chdir(ligGridDockDir)
    if not inducedFit:
        slurm.salloc("{}/glide {}/{}.in -WAIT".format(SCHRODINGER, ligGridDockDir, ligToGrid), 
            1, "10:00:00") #jobString, numProcessors, timeLimit
    else:
        slurm.salloc("{}/ifd -NGLIDECPU 1 -NPRIMECPU 1 {}/{}.inp -SUBHOST localhost -WAIT".format(SCHRODINGER, ligGridDockDir, ligToGrid),
                2, "50:00:00") #jobstring, numProcessors, timeLimit

    #Return (ligand,grid) so we know the directory to check for success
    return (ligand, grid)

def dockHelper(ligGrid):
    return dock(*ligGrid)

def dockDataset(dataset, xDock=True):
    global GLIDE_DIR
    if xDock:
        GLIDE_DIR = 'xglide'
    else:
        GLIDE_DIR = 'glide'

    gridsDir = "{}/{}/{}".format(DATA_DIR, dataset, GRIDS_DIR)

    structures = [o for o in os.listdir(gridsDir) if os.path.isdir(os.path.join(gridsDir,o))]

    dock("AR", "3B5R", "1XOW", False, True)

    toDock = [] #List of (ligand, grid) tuples that need to be submitted for docking

    for ligand in structures:
        for grid in structures:
            if not glideExists(dataset, ligand, grid):
                toDock.append((dataset, ligand, grid, xDock))

    print("{} ligand-grid pairs are missing! Checking for failures...".format(str(len(toDock))))
    
    dockFailures = [] #List of structures that failed because docking couldn't find a good pose

    for dataset, ligand, grid, xDock in toDock:
        if glideFailed(dataset, ligand, grid):
            dockFailures.append((dataset, ligand, grid, xDock))

    print("Of {} missing ligand-grid pairs, {} were failures (listed below). Not submitting these...".format(str(len(toDock)), str(len(dockFailures))))
    print([(x[1], x[2]) for x in dockFailures])

    toDock = [x for x in toDock if x not in dockFailures]
    print("Submitting:")
    print([(x[1], x[2]) for x in toDock])

    pool = Pool(5) #We have 5 Glide licenses

    while len(toDock) != 0:
        currentlyDocking = toDock
        toDock = []

        for finishedLigand, finishedGrid in pool.imap_unordered(dockHelper, currentlyDocking):
            #Check to make sure that glide worked successfully
            if not glideExists(dataset, finishedLigand, finishedGrid) and not glideFailed(dataset, finishedLigand, finishedGrid):
                toDock.append((dataset, finishedLigand, finishedGrid, xDock))

    #Now all the potentially successful jobs are finished, dock the failed structures

    #First, recalculate the failed jobs, more failures could have occured from our submissions
    dockFailures = []

    for ligand in structures:
        for grid in structures:
            if not glideExists(dataset, ligand, grid) and glideFailed(dataset, ligand, grid):
                dockFailures.append((dataset, ligand, grid, False, True)) #xDock = False, inducedFit = True
    
    print("Submitting the failed docking runs as induced fit docking jobs...")
    print(dockFailures)

    #Here, we don't want to attempt a resubmission, maybe include this later? Implementing a first pass
    pool.imap_unordered(dockHelper, dockFailures) 

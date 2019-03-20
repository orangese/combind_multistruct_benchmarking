import os
import sys
import glob
import numpy as np
import re

from containers import Protein

from schrodinger.structutils.analyze import evaluate_smarts_canvas
from schrodinger.structutils.rmsd import ConformerRmsd
from schrodinger.structure import StructureReader, StructureWriter

def pairwise_rmsds(poseviewer_path_1, poseviewer_path_2, rmsd_file, max_poses):
    """Calculates all pairwise RMSDs between all poses in two Glide poseviewer outputs

    Parameters
    ----------
    poseviewer_path_1 : str
        Path to poseviewer file 1, output by Glide docking
    poseviewer_path_2 : str
        Path to poseviewer file 2, output by Glide docking
    rmsd_file : str
        Path to file to write containing RMSDs between poses
    max_poses : int
        Number of top poses to use from each poseviewer file

    Returns
    -----------
    None
    """
    pv1 = list(StructureReader(poseviewer_path_1))[1:]
    pv2 = list(StructureReader(poseviewer_path_2))[1:]

    rmsds = {}
    for i, pose1 in enumerate(pv1[:max_poses]):
        for j, pose2 in enumerate(pv2[:max_poses]):
            calc = ConformerRmsd(pose1, pose2)
            calc.use_heavy_atom_graph = True
            rmsds[(i, j)] = calc.calculate()
        
    with open(rmsd_file, 'w') as f:
        for (i, j), rmsd in rmsds.items():
            f.write('{},{},{}\n'.format(i, j, rmsd))

def calc_scores(prot, ligand, structures):
    """Find the docking scores for a ligand to protein structures of interest
    
    Parameters
    ----------
    prot : str
        Name of protein (e.g. A2AR)
    ligand : str
        Name of ligand (e.g. NECA) - note: even though we name the ligands with
        "_lig" convention, we do not need to include "_lig" in this case
    structures : list of str
        List of strings of crystal structure names. We will return the docking energies
        of the ligand to these crystal structures. e.g. '["2YDO", "2YDO_181S"]'

    Returns
    -------
    allScores
        Dictionary of lists. Keys are structure names input in "strcutures" argument. 
        Values are list of docking scores.
    
    """

    ligand = ligand + '_lig'
    predict_data = Protein(prot) # prot is protein name
    allPoses = {}
    for st in structures:
        scores = []
        predict_data.load_docking([ligand], st = st)
        poses = predict_data.docking[st].ligands[ligand].poses
        for pose in poses:
            scores.append(pose.gscore)
        scores = np.array(scores)
        allPoses[st] = scores
    
    return(allPoses)
    
def process_rmsds(prot, ligToProcess):
    """Find all pairwise RMSDs between all docking runs that involve a ligand of interest for
       a given protein system.
       
    Parameters
    ----------
    prot : str
        Name of protein with docking results (e.g. A2AR)
    ligToProcess : str
        Name of ligand to get pairwise RMSDs for docking in all structures (e.g. NECA)
    
    Returns
    -------
    None
        But prints all pairwise RMSD comparison files into rmsd/ folder
    
    """

    resultsDir = '/scratch/PI/rondror/combind/mutagenesis/{}/docking/confgen_es4/'.format(prot)
    outputDir = '/scratch/PI/rondror/combind/mutagenesis/{}/rmsd/'.format(prot)
    subFolders = os.listdir(resultsDir)
    
    # choose only folders in resultsDir
    folders = [folder for folder in os.listdir(resultsDir) if os.path.isdir(os.path.join(resultsDir, folder))]
    
    # keep track of all ligands we have docked to
    ligs = set()
    
    for folder in folders:
        ligs.add(folder.split('_')[0])
    
    # rename folders to have full path
    folders = [os.path.join(resultsDir, folder) for folder in folders]
    
    pvFiles = []
    for lig in ligs:
        ligPvFiles = []
        for folder in folders:
            if '{}_lig'.format(lig) in folder:
                os.chdir(folder)
                pvFile = glob.glob('./*pv*')
                if not len(pvFile) > 0:    # if no poses found, skip checking for pvFile
                    os.chdir('../')
                    continue
                pvFile = pvFile[0]
                pvFile = os.path.join(folder, pvFile)
                ligPvFiles.append(pvFile)
                os.chdir('../')
        pvFiles.append(ligPvFiles)
    
    # though we gathered the poseviewer files for all ligands, we now only output RMSD comparison
    # for the ligand of interest. It may be worthwhile to change this so it outputs RMSD comparisons
    # for all ligands
    ligs = list(ligs)
    ligIdx = ligs.index(ligToProcess)
    ligPv = pvFiles[ligIdx]
    ligPv = sorted(ligPv, key=len)     # sort entries by length - i.e. prioritize comparisions relative
                                       # to native receptor first
    for i in range(len(ligPv)):
        for j in range(i+1, len(ligPv)):
            pv1 = ligPv[i]
            pv2 = ligPv[j]
            pv1Base = os.path.basename(pv1)[:-9] # drop _pv.maegz from name
            pv2Base = os.path.basename(pv2)[:-9] # drop _pv.maegz from name
            rmsdName = '{}-{}_rmsd.csv'.format(pv1Base, pv2Base)
            rmsdName = os.path.join(outputDir, rmsdName)
            if not os.path.exists(rmsdName):
                pairwise_rmsds(pv1, pv2, rmsdName, 100)
    
def all_pose_comparison(prot, refStruct, lig):
    """Collate all pose comparisons for a reference ligand + protein strcuture to the
    poses of the same reference ligand + other protein strcutures

    Parameters
    ----------
    prot : str
        Name of protein (e.g. A2AR)
    refStruct : str
        Which structure of the protein is the reference for the pairwise pose RMSD predictions?
        (e.g. 2YDO_181S)
    lig : str
        Which ligand are we collating pose pairwise predictions for? (e.g. NECA)

    Returns
    -------
    comparisonData : dictionary of lists
        Keys of dictionary are receptor names (e.g. 2YDO_182A), values are pairwise pose predictions
        between reference structure + ligand, and the structure + ligand
    
    """
    workingDir = '/scratch/PI/rondror/combind/mutagenesis/{}/rmsd'.format(prot)
    comparisonFiles = glob.glob(os.path.join(workingDir, '{}_lig-to-{}-*'.format(lig, refStruct)))
    comparisonData = {}
    for f in comparisonFiles:
        data = read_pairwise_poses_rmsd_file(f)
        receptorMatch = re.match(r'.*-([^-]*)_rmsd.csv$', f)
        if receptorMatch is None:
            raise Exception('File name is incorrect: {}'.format(f))
        else:
            receptor = receptorMatch.group(1)
        comparisonData[receptor] = data
    
    return(comparisonData)

def read_pairwise_poses_rmsd_file(fileName):
    """Read file containing RMSDs between poses from different structures
    
    Parameters
    ----------
    fileName : str
        Name of file containing RMSDs comparing between poses from different receptors
    
    Returns
    -------
    rmsdData : numpy array
    """

    ##### NOT YET tested
    
    import pdb; pdb.set_trace()
    pose1Max = -1
    pose2Max = -1
    with open(fileName, 'r') as f:
        for line in f:
            pose1, pose2, _ = line.rstrip().split(',')
            pose1Max = max(int(pose1), pose1Max)
            pose2Max = max(int(pose2), pose2Max)
            
    assert(pose1Max != -1)
    assert(pose2Max != -1)
    
    rmsdData = np.zeros((pose1Max+1, pose2Max+1))
    
    with open(fileName, 'r') as f:
        for line in f:
            pose1, pose2, rmsd = line.rstrip().split(',')
            rmsdData[int(pose1), int(pose2)] = float(rmsd)
        
    return(np.array(rmsdData))

def read_pose_to_native_rmsd_file(fileName):
    """Read file comparing docked poses for a given structure + ligand
    to native pose for the ligand
    
    Parameters
    ----------
    fileName : str
        Path to file containing RMSDs of docked poses to native pose
    
    Returns
    -------
    rmsdData : numpy array
    """
    rmsdData = []
    with open(fileName, 'r') as f:
        f.readline()
        for line in f:
            data = line.rstrip()
            data = data.split(',')
            data = [entry.strip('"') for entry in data]
            rmsdData.append(float(data[3]))
    return(np.array(rmsdData))

def average_pose_over_mutants(prot, lig, control=False):
    """Combine the docking results of a ligand across all mutants to generate improved predictions
    of the native docking conformation. First, find all poses of the ligand to the WT reference
    receptor. Next, find all poses of the ligand to all receptors. For each pose of ligand + WT ref,
    find the first pose of the ligand + other receptors within 2 Angstroms of the lig + WT ref pose. 
    Rescore the lig + WT ref pose with the average of the Glide scores for the first poses within 2
    Angstroms found among the other receptors. Output a re-ranked list of docking results.

    Parameters
    ----------
    lig : str
        Name of ligand to rescore (e.g. NECA)
    prot : str
        Name of protein of interest (e.g. A2AR)
    control : boolean, False by default
        Do we want to incorporate the mutant receptors that should have observed 
        binding (control=False)? Or do we want to incorporate the mutant receptors 
        that should have no observed binding (control=True)?

    Returns
    -------
    tuple
        (posesToNativeRmsd, rescoredPosesToNativeRmsd)
        posesToNativeRmsd : RMSD of best scored poses by Glide of ligand + WT receptor
        rescoredPosesToNativeRmsd : RMSD of poses by Glide of ligand + WT receptor, after rescoring
            and reordering the poses by incorporating docking to mutant receptors
    
    Will also output several files to the rescored/ folder:
        - file of the posesToNativeRmsd and rescoredPosesToNativeRmsd
        - file of the scores of each pose for all mutant receptor used in the analysis
    """

    refProt = ''
    if (prot == 'A2AR'):
        refProt = '2YDO'
    elif (prot == 'DHFR'):
        refProt = '3GHW'

    poseComparisons = all_pose_comparison(prot, refProt, lig)
    if (prot == 'A2AR'):
        if not control:
            structures = ['2YDO', '2YDO_181S', '2YDO_250F', '2YDO_250Y', '2YDO_277N', '2YDO_277T', '2YDO_281T']
        else:
            structures = ['2YDO', '2YDO_182A', '2YDO_250A', '2YDO_253A', '2YDO_278A', '2YDO_281A']
    elif (prot == 'DHFR'):
        structures = ['3GHW', '3GHW_22F', '3GHW_22W', '3GHW_22Y', '3GHW_35K', '3GHW_35K_64F', '3GHW_35S', '3GHW_64F', '3GHW_64S', '3GHW_35S_64S', '3GHW_64F', '3GHW_64S']
    allScores = calc_scores(prot, lig, structures)
    

    nativeRmsdDir = '/scratch/PI/rondror/combind/mutagenesis/{}/docking/confgen_es4/'.format(prot)
    posesToNativeRmsdsFile = os.path.join(nativeRmsdDir, '{}_lig-to-{}/rmsd.csv'.format(lig, refProt))
    posesToNativeRmsds = read_pose_to_native_rmsd_file(posesToNativeRmsdsFile)
    posesToNativeRmsds = np.array(posesToNativeRmsds)[:99]

    # find how pose 0 to pose 99 vary based on score
    # import pdb; pdb.set_trace()
    averagePoseScore = []
    referencePoseScores = allScores[refProt]
    mutantStructures = [struct for struct in structures if struct in poseComparisons.keys()]
    for p in np.arange(0, 99):
        poseScores = []
        poseScores.append(referencePoseScores[p])
        for mutant in mutantStructures:
            mutantScores = allScores[mutant]
            mutantPoseRmsds = poseComparisons[mutant]
            mutantPoseRmsds = mutantPoseRmsds[p] # look only at the pose of interest from the reference docking
            mutantPoseRmsds = mutantPoseRmsds < 2  # find all pose comparisons with RMSD < 2
            poseIdx = np.argmax(mutantPoseRmsds)
            poseScore = mutantScores[poseIdx]
            poseScores.append(poseScore)
        averagePoseScore.append(np.average(poseScores))

    poseOrder = np.argsort(averagePoseScore)
    reorderedPosesToNativeRmsds = np.array(posesToNativeRmsds[poseOrder])
    
    # write results to file
    fullResults = np.stack((posesToNativeRmsds, reorderedPosesToNativeRmsds), axis=1)
    outputDir = '/scratch/PI/rondror/combind/mutagenesis/{}/rescored'.format(prot)
    outputFile = lig + '_rescored'
    if control:
        outputFile += '_control.csv'
    else:
        outputFile += '.csv'
    np.savetxt(os.path.join(outputDir, outputFile), fullResults)
    
    # as a sanity check, store the average of the pose scores for the mutants
    allMutantScores = [allScores[mutant][:99] for mutant in mutantStructures]
    allMutantScores = np.stack(allMutantScores, axis=1)
    # allMutantScores = np.average(allMutantScores, axis=1)
    mutantScoreFile = lig + '_mutant_scores'
    if control:
        mutantScoreFile += '_control.csv'
    else:
        mutantScoreFile += '.csv'
    np.savetxt(os.path.join(outputDir, mutantScoreFile), allMutantScores)
    
    return(posesToNativeRmsds, reorderedPosesToNativeRmsds)


if __name__ == '__main__':
    # import pdb; pdb.set_trace()
    prog = sys.argv[1]
    if prog == 'rmsd':
        poseviewer_path_1, poseviewer_path_2, rmsd_file, max_poses = sys.argv[2:]
        pairwise_rmsds(poseviewer_path_1, poseviewer_path_2, rmsd_file, int(max_poses))
    elif prog == 'score':
        prot, ligand, structures = sys.argv[2:]
        structures = structures.split(',')
        calc_scores(prot, ligand, structures)
    elif prog == 'compare':
        prot, ligand = sys.argv[2:]
        process_rmsds(prot, ligand)
    elif prog == 'aggregate':
        prot, ligand, control = sys.argv[2:]
        if control.lower() == 'true':
            average_pose_over_mutants(prot, ligand, True)
        else:
            average_pose_over_mutants(prot, ligand, False)

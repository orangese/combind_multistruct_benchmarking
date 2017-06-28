import os
import sys
import tqdm
import itertools
from math import exp
import numpy as np
import os.path

os.chdir("../3_score/")
from fingerprint import FuzzyFingerPrint
from ligand import Ligand
from pose import Pose
from rmsd import RMSD


def load_data(data_set_dir, rmsd_file, ligands_dir, grids_dir, glide_dir, crystal_fp_file, docking_fp_dir):
    os.chdir(data_set_dir)

    # if not all structures (grids or ligands) are desired, remove them here
    gridstructs = [d.upper() for d in os.listdir(grids_dir) if os.path.isdir(os.path.join(grids_dir, d))]
    print 'Found ' + str(len(gridstructs)) + ' grids'
    ligstructs = [d.upper() for d in os.listdir(ligands_dir) if os.path.isfile(os.path.join(ligands_dir, d))]
    ligstructs = map(lambda x: x.split("_")[0], ligstructs)
    print 'Found ' + str(len(ligstructs)) + ' ligands'
    #Sort our ligands by the order in which we perform SAR analysis
    #all_gridstructs_upper = map(lambda x : x.upper(), gridstructs)
    ligstructs.sort(key=lambda lig: gridstructs.index(lig) if lig in gridstructs else len(gridstructs) + 1)
    
    crystals = load_crystals(crystal_fp_file, gridstructs)
    glides = load_glides(gridstructs, docking_fp_dir, glide_dir, rmsd_file)
    
    return (crystals, glides)


def load_crystals(crystal_fp_file, grids):
    print 'Loading crystal structures...'
    crystals = {} # PDB id : Pose
    for line in open(crystal_fp_file):
        if len(line) < 2:
            continue

        struct, ifp = line.strip().split(';')

        fp = FuzzyFingerPrint.compact_parser(ifp, struct.upper())
        crystals[struct.upper()] = Pose(0.0, fp, 0, 0)

    #Check to see if we are missing any crystal fingerprints
    crystalFPDiff = set( grids ).difference(set( crystals.keys() ))

    if(len(crystalFPDiff) != 0):
        print("Missing the following Crysal Fingerprints! Check " + crystal_fp_file)
        print(list(crystalFPDiff))
    return crystals

def load_gscores(gridstructs, glide_dir):
    print 'Loading glidescores...'
    # Ligand -> Grid Structure -> List of Glide Scores sorted by pose number (Lowest Scores to Highest Scores)
    gscores = {}

    failed_to_dock = 0
    total = 0

    for ligand in gridstructs:
        gscores[ligand] = {}
        for grid in gridstructs:
            total += 1
            gscores[ligand][grid] = []

            #Go through all possible permutations for file capitalization
            #NOTE: Files still have to be in the {}_ligand-to-{} format!
            cap_permutations = [(ligand, grid), (ligand.upper(), grid.upper()), (ligand.upper(), grid.lower()), (ligand.lower(), grid.upper()), (ligand.lower(), grid.lower())]
            s_file = None
            for tLigand, tGrid in cap_permutations:
                l_to_g = '{}_ligand-to-{}'.format(tLigand, tGrid)
                pv = glide_dir + l_to_g + '/' + l_to_g + '_pv.maegz'

                if os.path.exists(pv):
                    break
            else:
                #print 'ligand, grid', ligand, grid, 'did not dock.'
                failed_to_dock += 1
                continue

            s_file = open(glide_dir + l_to_g + '/' + l_to_g + '.rept')
            line = s_file.readline().strip().split()
            while not line or line[0] != '1' or (len(line) != 19 and (len(line) > 1 and line[1] != "1" and len(line) != 18)): # hack - why is title blank sometimes?
                line = s_file.readline().strip().split()
            while line:
                # Rank', 'Title', 'Lig#', 'Score', 'GScore'
                rank, title, lig, score, gscore = line[:5]
                gscores[ligand][grid] += [float(gscore)]
                line = s_file.readline().strip().split()
    if failed_to_dock == 0: print 'All ligands and structures successfully docked. Nice!'
    else: print str(failed_to_dock) + ' of ' + str(total) + ' ligand, structure pairs failed to dock.'
    return gscores


def load_rmsds(rmsd_file):
    print 'Loading rmsds...'
    # Ligand -> Grid Structure -> List of RMSDs ordered by GLIDE Pose Number
    rmsds = {}
    for line in open(rmsd_file):
        n, data = line.strip().split(':')
        ligand, grid = n.split('-to-')
        ligand = ligand.split('_')[0].upper()
        grid = grid.split('_')[0].upper()

        if grid not in rmsds: rmsds[grid] = {}
        rmsds[grid][ligand] = RMSD.read(data)
    return rmsds


def load_glides(gridstructs, docking_fp_dir, glide_dir, rmsd_file):
    print 'Loading docking results...'
    gscores = load_gscores(gridstructs, glide_dir)
    rmsds = load_rmsds(rmsd_file)
            
    glides = {}
    for ligand in gridstructs:
        glides[ligand] = {}
        for grid in gridstructs:
            glides[ligand][grid] = Ligand(None)
            if len(gscores[ligand][grid]) == 0:
                # put in one filler pose
                glides[ligand][grid].add_pose(Pose(100,None,0,0),0)
                continue # ligand, grid, did not dock.

            #Go through all possible permutations for file capitalization
            #NOTE: Files still have to be in the {}_ligand-to-{} format!            
            cap_permutations = [(ligand, grid), (ligand.upper(), grid.upper()), 
                                (ligand.upper(), grid.lower()), (ligand.lower(), grid.upper()), (ligand.lower(), grid.lower())]
            for tLigand, tGrid in cap_permutations:
                fnm = docking_fp_dir + "{}_ligand-to-{}.fp".format(tLigand, tGrid, tLigand, tGrid)
                try:
                    s_file = open(fnm)
                    break
                except:
                    pass

            for pose_num, line in enumerate(s_file):
                try:
                    ifp = line.strip()
                except Exception as e:
                    print("Error: No fingerprint generated for pose " + line.strip() + " in " + fnm)
                    break

                #if pose_num < 50: #ONLY IMPORT 50 POSES, REMOVE LATER
                try:
                    glides[ligand][grid].add_pose(Pose(rmsds[grid][ligand].get_rmsd(pose_num),
                                                       FuzzyFingerPrint.compact_parser(ifp, ligand), pose_num,
                                                       gscores[ligand][grid][pose_num]), pose_num)
                except Exception as e:
                    print 'key error -- check capitalization of grids, ligands', e
                    print rmsds.keys()
                    print gscores.keys()
                    print rmsds[grid].keys()
                    print gscores[ligand].keys()
                    break
    return glides


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
#from cluster import Cluster
from pose import Pose
from rmsd import RMSD


def load_data(data_set_dir, rmsd_file, ligands_dir, grids_dir, glide_dir, crystal_fp_file, docking_fp_dir):
    os.chdir(data_set_dir)

    # if not all structures (grids or ligands) are desired, remove them here
    gridstructs = [d for d in os.listdir(grids_dir) if os.path.isdir(os.path.join(grids_dir, d))]

    ligstructs = [d for d in os.listdir(ligands_dir) if os.path.isfile(os.path.join(ligands_dir, d))]
    ligstructs = map(lambda x: x.split("_")[0], ligstructs)

    #Sort our ligands by the order in which we perform SAR analysis
    all_gridstructs_upper = map(lambda x : x.upper(), gridstructs)
    ligstructs.sort(key=lambda lig: all_gridstructs_upper.index(lig.upper()) if lig.upper() in all_gridstructs_upper else len(all_gridstructs_upper) + 1)
    
    crystals = load_crystals(crystal_fp_file)#, all_gridstructs_upper)
    glides = load_glides(gridstructs, docking_fp_dir, glide_dir, rmsd_file)
    
    return (crystals, glides)


def load_crystals(crystal_fp_file):#, all_gridstructs_upper):
    crystals = {} # PDB id : Pose
    for line in open(crystal_fp_file):
        if len(line) < 2:
            continue

        struct, ifp = line.strip().split(';')

        #if struct.upper() in all_gridstructs_upper:
        fp = FuzzyFingerPrint.compact_parser(ifp, struct.upper())
        crystals[struct.upper()] = Pose(0.0, fp, 0, 0)

    #Check to see if we are missing any crystal fingerprints
    #crystalFPDiff = set(all_gridstructs_upper).difference(set(map(lambda x: x.upper(), crystals.keys())))

    #if(len(crystalFPDiff) != 0):
    #    print("Missing the following Crysal Fingerprints! Check " + crystal_fp_file)
        #print(list(crystalFPDiff))
        #exit()
    return crystals

def load_gscores(gridstructs, glide_dir):
    # Ligand -> Grid Structure -> List of Glide Scores sorted by pose number (Lowest Scores to Highest Scores)
    gscores = {}

    for ligand in gridstructs:
        gscores[ligand] = {}
        for grid in gridstructs:
            gscores[ligand][grid] = []

            #Go through all possible permutations for file capitalization
            #NOTE: Files still have to be in the {}_ligand-to-{} format!
            cap_permutations = [(ligand, grid), (ligand.upper(), grid.upper()), (ligand.upper(), grid.lower()), (ligand.lower(), grid.upper()), (ligand.lower(), grid.lower())]
            s_file = None
            for tLigand, tGrid in cap_permutations:
                fnm = glide_dir + "{}_ligand-to-{}/{}_ligand-to-{}.rept".format(tLigand, tGrid, tLigand, tGrid)
                try:
                    s_file = open(fnm)
                    break
                except:
                    pass
            if s_file == None:
                print 'Did not find ', tLigand, tGrid
                continue
            line = s_file.readline().strip().split()
            while not line or line[0] != '1' or (len(line) != 19 and (len(line) > 1 and line[1] != "1" and len(line) != 18)): # hack - why is title blank sometimes?
                line = s_file.readline().strip().split()
            while line:
                # Rank', 'Title', 'Lig#', 'Score', 'GScore'
                rank, title, lig, score, gscore = line[:5]
                gscores[ligand][grid] += [float(gscore)]
                line = s_file.readline().strip().split()
    return gscores


def load_rmsds(rmsd_file):
    # Ligand -> Grid Structure -> List of RMSDs ordered by GLIDE Pose Number
    rmsds = {}
    #for line in open('rmsd_table.csv'):
    for line in open(rmsd_file):
        n, data = line.strip().split(':')
        ligand, grid = n.split('-to-')
        ligand = ligand.split('_')[0].upper()
        grid = grid.split('_')[0].upper()

        if grid not in rmsds: rmsds[grid] = {}
        rmsds[grid][ligand] = RMSD.read(data)
    return rmsds


def load_glides(gridstructs, docking_fp_dir, glide_dir, rmsd_file):
    gscores = load_gscores(gridstructs, glide_dir)
    rmsds = load_rmsds(rmsd_file)
            
    glides = {}
    for ligand in gridstructs:
        glides[ligand] = {}
        for grid in gridstructs:
            glides[ligand][grid] = Ligand(None)

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
                if pose_num > 50:
                    break
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


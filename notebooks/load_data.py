import os
import sys
import os.path

os.chdir("../3_score/")
from fingerprint import FuzzyFingerPrint
from ligand import Ligand
from pose import Pose
from rmsd import RMSD


def load_ligand_info(data_set_dir):
    os.chdir(data_set_dir)
    
    lig_info = {}
    for line in open(data_set_dir + 'lig_sizes.txt'):
        a = line.strip().split(',')
        lig_info[a[0].split('/')[-1][:4]] = [float(i) for i in a[1:]]
    return lig_info

def load_data(data_set_dir, rmsd_file, glide_dir, crystal_fp_file, docking_fp_dir, w=[10,10,10,0,1]):
    os.chdir(data_set_dir)

    crystals = load_crystals(crystal_fp_file, w=w)

    gscores, rmsds = load_gscores(glide_dir, rmsd_file)
    ifps = load_ifps(docking_fp_dir, w=w)
       
    glides = load_glides(gscores, rmsds, ifps)
    
    return (crystals, glides)


def load_crystals(crystal_fp_file, w=None):
    print 'Loading crystal structures...'
    crystals = {} # PDB id : Pose
    
    if not os.path.isfile(crystal_fp_file):
        print 'Crystal structure fingerprint file not found.'
    else:
        fp_file = open(crystal_fp_file)
        for line in fp_file:
            struct, ifp = line.strip().split(';')
            crystals[struct] = Pose(0.0, FuzzyFingerPrint.compact_parser(ifp, struct, w=w), 0, 0)
        fp_file.close()
    return crystals

def load_gscores(glide_dir, rmsd_file):
    print 'Loading glidescores...'
    # Ligand -> Grid Structure -> List of Glide Scores sorted by pose number (Lowest Scores to Highest Scores)
    (total, undocked) = 0, 0
    gscores = {}
    rmsds = {}
    for gdir in os.listdir(glide_dir):
        if '_ligand-to-' not in gdir: continue
        lig, struct = gdir.split('_ligand-to-')
        total += 1
        if not (os.path.exists(glide_dir + gdir + '/' + gdir + '_pv.maegz') and os.path.exists(glide_dir + gdir + '/' + gdir + '.rept')):
            undocked += 1
            continue

        if lig not in gscores: gscores[lig] = {}
        gscores[lig][struct] = []

        gscore_file = open(glide_dir + gdir + '/' + gdir + '.rept')
        line = gscore_file.readline().strip().split()
        while not line or line[0] != '====':
            line = gscore_file.readline().strip().split()
        line = gscore_file.readline().strip().split()
        
        while line:
            # Rank', 'Title', 'Lig#', 'Score', 'GScore'
            if line[2] == '1': gscores[lig][struct].append(float(line[4]))
            elif line[1] == '1': gscores[lig][struct].append(float(line[3]))
            else:
                print 'Lig# 1 not found. quitting. lig, grid: ', lig, struct
                break
            #print line
            if line[-1] != '--':
                if struct not in rmsds: rmsds[struct] = {}
                if lig not in rmsds[struct]: rmsds[struct][lig] = []
                rmsds[struct][lig].append(float(line[-1]))
            #    print line[-1]
            line = gscore_file.readline().strip().split()
        gscore_file.close()
    if len(rmsds.keys()) == 0:
        rmsds = load_rmsds(rmsd_file)
    print '{} of {} total pairs failed to dock.'.format(undocked, total)
    return gscores, rmsds

def load_rmsds(rmsd_file):
    print 'Loading rmsds...'
    # Ligand -> Grid Structure -> List of RMSDs ordered by GLIDE Pose Number
    rmsds = {}
    all_rmsds = open(rmsd_file)
    for line in all_rmsds:
        n, data = line.strip().split(':')
        lig, struct = n.split('_ligand-to-')

        if struct not in rmsds: rmsds[struct] = {}
        rmsds[struct][lig] = RMSD.read(data).data
    all_rmsds.close()
    return rmsds

def load_ifps(docking_fp_dir, w=None):
    print 'Loading fingerprints...'
    ifps = {}
    if not os.path.exists(docking_fp_dir):
        print 'Glide fingerprint directory not found.'
    else:
        for fnm in [i for i in os.listdir(docking_fp_dir) if i[-3:] == '.fp']:
            lig, struct = fnm[:-3].split('_ligand-to-')

            if struct not in ifps: ifps[struct] = {}
            ifps[struct][lig] = []

            fp_file = open(docking_fp_dir + fnm)
            for pose_num, line in enumerate(fp_file):
                ifps[struct][lig].append(FuzzyFingerPrint.compact_parser(line.strip(), lig, w=w))
            fp_file.close()
    return ifps

def load_glides(gscores, rmsds, ifps):
    print 'Loading docking results...'
    glides = {}
    for lig in gscores.keys():
        glides[lig] = {}
        for struct in gscores[lig].keys():
            if len(gscores[lig][struct]) == 0:
                continue # ligand, grid, did not dock.    

            try:
                assert len(gscores[lig][struct]) == len(rmsds[struct][lig])
                #assert len(ifps[struct][lig]) <= len(rmsds[struct][lig])
            except:
                print 'lig, struct: ', lig, struct, rmsds.keys()
                print rmsds[struct].keys()
                print len(rmsds[struct][lig])
                print 'gscores, rmsds, ifps: ', len(gscores[lig][struct]), len(rmsds[struct][lig]), len(ifps[struct][lig])

            glides[lig][struct] = Ligand(None)
            for pose in range(len(gscores[lig][struct])):
                rmsd = rmsds[struct][lig][pose]
                gscore = gscores[lig][struct][pose]
                ifp = ifps[struct][lig][pose] if ifps != {} and pose < len(ifps[struct][lig]) else None
                glides[lig][struct].add_pose(Pose(rmsd, ifp, pose, gscore), pose)

    return glides


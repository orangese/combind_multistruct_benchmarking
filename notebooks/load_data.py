import os
import sys
import os.path

os.chdir("../3_score/")
from fingerprint import FuzzyFingerPrint
from ligand import Ligand
from pose import Pose
from rmsd import RMSD

MAX_NUM_POSES = 500

# grouped in order of biggest ligands to smallest
similarity_groups = {
    'B1AR_all': {
        0: ['2Y00','2Y01','2Y02','4AMI','4AMJ'],
        1: ['2VT4','2YCW','2YCX','2YCY','2YCZ','4BVN','5F8U'],
        2: ['2Y03','2Y04'], # these are substructures of group 0
        3: ['3ZPQ','3ZPR'], # these are only slightly smaller than 2
        4: ['4GPO']         # no ligand bound
    },
    'B2AR_all': {
        0: ['4LDL','4QKX','4LDE','3SN6','3PDS','3P0G'],
        1: ['5JQH','5D6L','5D5B','5D5A','4GBR','2RH1','3NYA','3NY8','3NY9'],
        2: ['4LDO','3D4S'],
        3: ['2R4S','2R4R','3KJ6']
    }
}

def sort_data(receptor, crystals, glides, combine_structs=False, exclude_duplicates=False):
    if glides == {}: return crystals, {}, crystals.keys(), crystals.keys()
    exclude = {
        'B1AR_all':['2Y01','4BVN','5F8U'],
        'B2AR':['3P0G'],
        'HSP90':['4YKY','4YKZ','4YKX']
    }
    
    ligs = sorted([i for i in glides.keys()])
    structs = []
    for l in ligs:
        structs.extend([i for i in glides[l] if i not in structs])# and i not in exclude.get(receptor,[])])
    structs.sort()
    if receptor in similarity_groups:
        sg = similarity_groups[receptor]
        f = lambda x: [i for i in range(max(sg.keys())+1) if x in sg[i]][0]
        ligs.sort(key=f)
        structs.sort(key=f)
    
    if combine_structs:
        for l in ligs:
            all_poses = []
            for s in glides[l]:
                all_poses.extend([pose for num, pose in glides[l][s].poses.items()])
            all_poses.sort(key=lambda x: x.gscore)
            glides[l]['all'] = Ligand(None)
            glides[l]['all'].poses = {i: p for i, p in enumerate(all_poses[:MAX_NUM_POSES])}
        structs.append('all')

    new_glides = {l:{s:glides[l][s] for s in glides[l] if s in structs} for l in glides if l in ligs}
    new_crystals = {l:crystals[l] for l in crystals if l in ligs}

    return new_crystals, new_glides, ligs, structs

def load_data(receptor,
        rmsd_file = 'xrmsd.csv',
        glide = 'xglide',
        crystal_ifp = 'ifp/xcrystal_ifp_3',
        glide_ifp = 'ifp/xglide_ifp_3',
        w = [10,10,10,1,0],
        require_fp=True,
        combine_structs=False,
        load_docking=True):

    data_set_dir = '/scratch/PI/rondror/docking_data/'+receptor
    glide_dir = '{}/{}/'.format(data_set_dir, glide)
    crystal_fp_dir = '{}/{}/'.format(data_set_dir, crystal_ifp)
    glide_fp_dir = '{}/{}/'.format(data_set_dir, glide_ifp)

    os.chdir(data_set_dir)
    
    crystals = load_crystals(crystal_fp_dir, w=w)

    glides = {}
    if load_docking:
        gscores, rmsds = load_gscores(glide_dir, rmsd_file)
        ifps = load_ifps(glide_fp_dir, w=w)
       
        glides = load_glides(gscores, rmsds, ifps, require_fp)
    
    return sort_data(receptor, crystals, glides, combine_structs)


def load_crystals(crystal_fp_dir, w=None):
    print 'Loading crystal structures...'
    crystals = {} # PDB id : Pose
    
    for f in [i for i in os.listdir(crystal_fp_dir) if i.split('.')[-1] == 'fp']:
        with open('{}/{}'.format(crystal_fp_dir, f)) as fp_file:
            for line in fp_file:
                struct, ifp = line.strip().split(';')
                if ifp == '':
                    continue
                crystals[struct] = Pose(0.0, FuzzyFingerPrint.compact_parser(ifp, struct, w=w), 0, 0, struct, struct)
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
        
        while line and len(gscores[lig][struct]) < MAX_NUM_POSES:
            # Rank', 'Title', 'Lig#', 'Score', 'GScore'
            if line[2] == '1': gscores[lig][struct].append(float(line[4]))
            elif line[1] == '1': gscores[lig][struct].append(float(line[3]))
            else:
                print 'Lig# 1 not found. quitting. lig, grid: ', lig, struct
                break
            if line[-1] != '--':
                if struct not in rmsds: rmsds[struct] = {}
                if lig not in rmsds[struct]: rmsds[struct][lig] = []
                rmsds[struct][lig].append(float(line[-1]))
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
                if pose_num >= MAX_NUM_POSES: break
                if line == '': continue
                ifps[struct][lig].append(FuzzyFingerPrint.compact_parser(line.strip(), lig, w=w))
            fp_file.close()
    return ifps

def load_glides(gscores, rmsds, ifps, require_fp):
    print 'Loading docking results...'
    glides = {}
    for lig in gscores.keys():
        glides[lig] = {}
        for struct in gscores[lig].keys():
            if require_fp and (lig not in ifps or struct not in ifps[lig] or len(ifps[lig][struct]) == 0):
                continue # ligand, grid, did not fingerprint
            elif not require_fp and len(gscores[lig][struct]) == 0:
                continue # did not dock

            glides[lig][struct] = Ligand(None)
            for pose in range(len(gscores[lig][struct])):
                rmsd = rmsds[struct][lig][pose]
                gscore = gscores[lig][struct][pose]
                try: ifp = ifps[struct][lig][pose]
                except: ifp = None
                glides[lig][struct].add_pose(Pose(rmsd, ifp, pose, gscore, struct, lig), pose)

    return glides


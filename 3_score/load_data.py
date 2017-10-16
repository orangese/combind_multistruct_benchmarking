import os
import sys
import os.path
import numpy as np

from pose import Pose

def load_residue_maps(receptor):
    r_map = {}
    os.chdir('/scratch/PI/rondror/docking_data/{}/rmsds'.format(receptor))
    for f in os.listdir('.'):
        assert f.split('.')[-1] == 'csv', f
        s1, s2 = f.split('.')[0].split('-')
        assert s2 not in r_map
        r_map[s2] = {}
        if s1 not in r_map:
            r_map[s1] = {}
        with open(f, 'r') as rmsd_f:
            for line in rmsd_f:
                if line.strip() == '' or line.split(',')[0] == 'Residue 1': continue
                try:
                    r1, r2, rmsd, ca_rmsd = line.strip().split(',')
                except:
                    print receptor, f, line
                    raise Exception()
                r_map[s1][r1] = r1
                r_map[s2][r2] = r1
    return r_map

def load_grids():
    grids = {}
    with open('{}/pdbbind_first_grids.txt'.format(data)) as f:
        for line in f:
            r, g = line.strip().split()
            grids['pdbbind_combo/'+r] = [g]

    with open('{}/pdbbind_second_grids.txt'.format(data)) as f:
        for line in f:
            r, g = line.strip().split()
            grids['pdbbind_combo/'+r].append(g)
    return grids

def load_2d_fp(receptor, crystals):
    fp2d = {}
    all_features = set([])
    for lig in os.listdir('/scratch/PI/rondror/docking_data/{}/ifp/2d_fp'.format(receptor)):
        with open('/scratch/PI/rondror/docking_data/{}/ifp/2d_fp/{}'.format(receptor, lig)) as f:
            for line in f:
                fp2d[lig.split('_')[0]] = line.strip().split(',')
                break
        for n in fp2d[lig.split('_')[0]]:
            if n not in all_features:
                all_features.add(n)

    f_list = list(all_features)
    return {l:[1 if f in fp2d[l] else 0 for f in f_list] for l in fp2d if l in crystals}

def load_combinations(receptor, structure, output_dir='scoring_output2'):
    all_scores = {}
    out_dir = '/scratch/PI/rondror/docking_data/{}/{}/'.format(receptor, output_dir)
    if not os.path.exists(out_dir):
        #print out_dir, 'does not exist'
        out_dir = '/scratch/PI/rondror/docking_data/{}/scoring_output'.format(receptor)

    for out_f in os.listdir(out_dir):
        s, k = out_f.split('.')[0].split('_struct_')
        if s != structure: continue
        k = int(k)
        #if s not in all_scores: all_scores[s] = {}
        if k not in all_scores: all_scores[k] = {}
        with open('{}/{}'.format(out_dir, out_f), 'r') as f:
            lag = 0
            ligs = None
            for line in f:
                if line[0] == '#': 
                    lag = 0
                    continue
                if lag % 2 == 0:
                    ligs = tuple(line.strip().split(','))
                else:
                    assert ligs is not None
                    poses = tuple([int(i) for i in line.strip().split(',')])
                    assert len(poses) == len(ligs)
                    all_scores[k][ligs] = poses
                lag += 1
    return all_scores

def load_data(receptor,
              struct,
              glide = 'xglide',
              crystal_ifp = 'xcrystal8',
              glide_ifp = 'xglide8',
              w = [0,0,1,1,1,0.5,1,1,1,0.02,0,0,0]):

    data_set_dir = '/scratch/PI/rondror/docking_data/'+receptor
    glide_dir = '{}/{}'.format(data_set_dir, glide)
    crystal_fp_dir = '{}/ifp/{}/'.format(data_set_dir, crystal_ifp)
    glide_fp_dir = '{}/ifp/{}/'.format(data_set_dir, glide_ifp)
    
    os.chdir(data_set_dir)
    #print receptor, data_set_dir, os.listdir('.')
    all_ligs = [i.split('_')[0] for i in os.listdir('{}/final_ligands'.format(data_set_dir))]
    
    r_map = load_residue_maps(receptor)

    if struct == 'crystals_only': 
        return load_crystals(all_ligs, crystal_fp_dir, w, r_map)

    gscores, rmsds = load_gscores(all_ligs, glide_dir, struct)
    ifps = load_ifps(glide_fp_dir, gscores.keys(), struct, w, r_map)
       
    return load_glides(gscores.keys(), struct, gscores, rmsds, ifps)

def parse_fp(fp, w, struct=None, r_map=None):
    fp_dict = {}
    for i in fp.split(';'):
        fp_prod = np.multiply(w, [float(j) for j in i.split(',')[1:]])
        if any(fp_prod):
            #fp_prod[0] += fp_prod[2]
            #fp_prod[1] += fp_prod[3]
            #fp_prod[2] = 0
            #fp_prod[3] = 0
            #fp_prod[9] = w[9]*(fp_prod[9]/float(w[9]))**0.5
            fp_dict[i.split(',')[0]] = fp_prod
#    if r_map: assert struct in r_map, '{} not found in rmap'.format(struct)
    if r_map and struct in r_map: 
        fp_dict = {r_map[struct].get(r, r):fp_dict[r] for r in fp_dict} 
    for r in fp_dict:
        if '(Zn)' in r:
            fp_dict['(Zn)'] = fp_dict[r]
            del fp_dict[r]
    return fp_dict

def load_crystals(all_ligs, crystal_fp_dir, w, r_map):
    crystals = {} # PDB id : fp
    for f in all_ligs:
        assert os.path.exists('{}/{}.fp'.format(crystal_fp_dir, f)), '{}: {} fp not found'.format(f, crystal_fp_dir)
        with open('{}/{}.fp'.format(crystal_fp_dir, f)) as fp_file:
            for line in fp_file:
                if line.strip() == '':
                    continue
                crystals[f] = parse_fp(line.strip(), w, f, r_map)
    return crystals

def load_gscores(all_ligs, glide_dir, struct):
    gscores = {lig:[] for lig in all_ligs}
    rmsds = {lig:[] for lig in all_ligs}
    for lig in all_ligs:
        gdir = '{}_ligand-to-{}'.format(lig, struct)
        #assert os.path.exists('{}/{}/{}_pv.maegz'.format(glide_dir, gdir, gdir)), '{}/{}/{} not docked'.format(glide_dir, gdir, gdir)
        #if not os.path
        if not os.path.exists('{}/{}/{}_pv.maegz'.format(glide_dir, gdir, gdir)):
            del gscores[lig]
            del rmsds[lig]
            continue
        with open('{}/{}/{}.rept'.format(glide_dir, gdir, gdir)) as gscore_file:
            prev = 0
            for line in gscore_file:
                line = line.strip().split()
                if len(line) < 1 or line[0] != str(prev + 1) or (line[1] != lig and line[1] != '{}_ligand'.format(lig)): continue 
                prev += 1
                if line[2] == '1': gscores[lig].append(float(line[4]))
                elif line[1] == '1': gscores[lig].append(float(line[3]))
                if line[-1] != '--': rmsds[lig].append(float(line[-1]))
    return gscores, rmsds

def load_ifps(docking_fp_dir, all_ligs, struct, w, r_map):
    ifps = {lig:[] for lig in all_ligs}
    for lig in all_ligs:
        fnm = '{}_ligand-to-{}.fp'.format(lig, struct)
        assert os.path.exists('{}/{}'.format(docking_fp_dir, fnm)), '{}/{}: no fp'.format(docking_fp_dir, fnm)

        with open('{}/{}'.format(docking_fp_dir, fnm)) as fp_file:
            for pose_num, line in enumerate(fp_file):
                if line == '': continue
                ifps[lig].append(parse_fp(line.strip(), w, struct, r_map))
    return ifps

def load_glides(all_ligs, struct, gscores, rmsds, ifps):
    glides = {lig:[] for lig in all_ligs}
    for lig in all_ligs:
        if len(gscores[lig]) == 0:
            del glides[lig]
            continue
        for pose in range(len(gscores[lig])):
            try: ifp = ifps[lig][pose]
            except: ifp = None
            glides[lig].append(Pose(rmsds[lig][pose], ifp, pose, gscores[lig][pose], struct, lig))
    return glides


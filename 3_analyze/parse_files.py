import os

def parse_fp_file(fp_file):
    ifps = {}
    try:
        with open(fp_file) as f:
            pose_num = 0
            for line in f:
                if line.strip() == '': continue
                if line[:4] == 'Pose':
                    pose_num = int(line.strip().split(' ')[1])
                    ifps[pose_num] = {}
                    continue
                sc_key, sc = line.strip().split('=')
                i,r,ss = sc_key.split('-')
                i = int(i)
                sc = float(sc)
                prev_sc = ifps[(i, r)] if (i,r) in ifps[pose_num] else 0
                ifps[pose_num][(i,r)] = max(prev_sc, sc)

    except Exception as e:
        print(e)
        print(fp_file, 'fp not found')
    if len(ifps) == 0:
        print('check', fp_file)
        return {}
    return ifps

def parse_glide_output(glide_dir):
    if not os.path.exists(glide_dir):
        return [], [], []
    pair = glide_dir.split('/')[-1]
    if os.path.exists('{}/{}_pv.maegz'.format(glide_dir, pair)):
        return parse_rept_file(glide_dir, pair)
    else:
        print('not finished', glide_dir)
        return [], [], []

def parse_rept_file(glide_dir, pair):
    lig, prot = pair.split('-to-')
    rept_file = '{}/{}.rept'.format(glide_dir, pair)
    rmsd_file = '{}/rmsd.csv'.format(glide_dir, pair)

    gscores, emodels, rmsds = [], [], []
    with open(rept_file) as fp:
        for line in fp:
            line = line.strip().split()
            if len(line) <= 1 or (line[1] != lig and line[1] != lig+'_out' and line[1] != '1'): continue
            rank, lig_name, lig_index, score = line[:4]
            emodel = line[13]
            if line[1] == '1':
                rank, lig_index, score = line[:3]
                emodel = line[12]
            gscores.append(float(score))
            emodels.append(float(emodel))
            
    if not os.path.exists(rmsd_file):
        return gscores, emodels, [None]*len(gscores)
    with open(rmsd_file) as fp:
        for line in fp:
            line = line.strip().split(',')
            if line[3] == '"RMSD"': continue
            rmsds.append(float(line[3][1:-1]))

    return gscores, emodels, rmsds

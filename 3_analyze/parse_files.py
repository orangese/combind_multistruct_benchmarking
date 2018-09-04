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

def parse_glide_output(g_dir):
    if not os.path.exists(g_dir):
        #print 'did not dock', g_dir
        return [], []
    os.chdir(g_dir)
    pair = g_dir.split('/')[-1]
    if os.path.exists('{}-out.maegz'.format(pair)):
        return parse_if_file(pair)
    elif os.path.exists('{}_pv.maegz'.format(pair)):
        return parse_rept_file(pair)
    else:
        print('not finished', g_dir)
        return [], []

def parse_if_file(pair):
    gscores = []
    rmsds = []
    try:
        with open('{}.rmsd'.format(pair)) as rmsd_f:
            for line in rmsd_f:
                pnum, rmsd = line.strip().split(' ')
                rmsds.append(float(rmsd))
        with open('{}_workdir/scoring_dir/report.csv'.format(pair)) as g_f:
            for i, line in enumerate(g_f):
                if i == 0: continue
                gscores.append(float(line.strip().split(',')[2])) # ifd score (protein and ligand)
    except:
        return [], []
    return gscores, rmsds

def parse_rept_file(pair):
    gscores, emodels, rmsds = [], [], []
    lig, prot = pair.split('-to-')
    with open('{}.rept'.format(pair)) as f:
        for line in f:
            line = line.strip().split()
            if len(line) <= 1 or (line[1] != lig and line[1] != lig+'_out' and line[1] != '1'): continue
            rank, lig_name, lig_index, score = line[:4]
            emodel = line[13]
            if line[1] == '1':
                rank, lig_index, score = line[:3]
                emodel = line[12]
            gscores.append(float(score))
            emodels.append(float(emodel))
            
    if not os.path.exists('rmsd.csv'):
        return gscores, emodels, [None]*len(gscores)
    with open('rmsd.csv') as f:
        for line in f:
            line = line.strip().split(',')
            if line[3] == '"RMSD"': continue
            #print line
            rmsds.append(float(line[3][1:-1]))

    return gscores, emodels, rmsds


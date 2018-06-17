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
                #if ss == '':
                #    ss = -1
                #else:
                #    ss = min([int(k) for k in ss.split(',')]) # min is the more conserved ss
                #if r_map and struct in r_map:
                #    r = r_map[struct].get(r, r)

                #ifps[pose_num][(i,r,'0')] = ifps[pose_num].get((i,r,'0'),0) + float(sc)*w[i]
                #if ss != '0':
                if (i,r) not in ifps[pose_num]:
                    ifps[pose_num][(i,r)] = float(sc)#{}
                #ifps[pose_num][(i,r)][ss] = float(sc)*w[i]

    except Exception as e:
        print e#, ss
        print fp_file, 'fp not found'
        #ifps.append(None)
    if len(ifps) == 0:
        print 'check', fp_file
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
        print 'not finished', g_dir
        return [], []
    #assert len(gscores) == len(rmsds)
    #return [Pose(rmsds[i], gscores[i]) for i in range(len(gscores))]

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
    gscores = []
    rmsds = []
    lig, prot = pair.split('-to-')
    with open('{}.rept'.format(pair)) as f:
        for line in f:
            line = line.strip().split()
            if len(line) <= 1 or (line[1] != lig and line[1] != lig+'_out' and line[1] != '1'): continue
            #if line[1] == lig:
            rank, lig_name, lig_index, score = line[:4]
            if line[1] == '1':
                rank, lig_index, score = line[:3]
            #rmsd = line[-1]

            gscores.append(float(score))

            #if rmsd == '--': rmsd = None
            #else: rmsd = float(rmsd)
            #rmsds.append(rmsd)
    if not os.path.exists('rmsd.csv'):
        return gscores, [None]*len(gscores)
    with open('rmsd.csv') as f:
        for line in f:
            line = line.strip().split(',')
            if line[3] == '"RMSD"': continue
            #print line
            rmsds.append(float(line[3][1:-1]))

    return gscores, rmsds


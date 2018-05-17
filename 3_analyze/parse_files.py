import os

def parse_mcss_size(mcss_path, ligs1, ligs2, struct):
    mcss_sizes = {}
    #print mcss_path
    if not os.path.exists(mcss_path): return {}
    for fname in sorted(os.listdir(mcss_path)):
        #print fname
        if fname.split('-to-')[-1] != '{}.csv'.format(struct): continue# or fname[0] == '.': continue
        l1, l2 = fname.split('-to-')[0].split('-')
        #print fname, l1, l2
        if not (l1 in ligs1 and l2 in ligs2):
            if not (l1 in ligs2 and l2 in ligs1):
                continue
        size = None
        try:
            with open('{}/{}'.format(mcss_path, fname)) as f:
                for i, line in enumerate(f):
                    line = line.strip().split(',')
                    if i == 0:
                        s1, s2, s3 = [0 if i == 'None' else float(i) for i in line]
                        mcss_sizes[(l1,l2)] = (s1,s2,s3)
                        break
        except:
            print 'mcss size error', fname

    return mcss_sizes

def parse_mcss(mcss_path, num_poses, all_pairs, min_size=10):
    mcss_scores = {}
    if not os.path.exists(mcss_path): return {}
    for l1,l2,st in all_pairs:
        fpath = '{}/{}-{}-to-{}.csv'.format(mcss_path, l1, l2, st)
        if not os.path.exists(fpath): continue
        size = None
        try:
            with open(fpath) as f:
                for i, line in enumerate(f):
                    line = line.strip().split(',')
                    if i == 0:
                        if line[-1] == 'None': break
                        s1, s2, s3 = [float(i) for i in line]
                        if int(s3) < min_size: break
                        size = s3/(s1+s2-s3)
                        mcss_scores[(l1, l2)] = {}
                    p1,p2,rmsd = int(line[0]), int(line[1]), float(line[2])
                    if rmsd*size > 100:
                        print l1, l2, p1, p2, rmsd, size
                    else:
                        mcss_scores[(l1, l2)][(p1, p2)] = rmsd#*size
                else:
                    if p1+1 < min(num_poses[l1], 100) or p2+1 < min(num_poses[l2], 100) or rmsd == float(10000):
                        print 'hmmmm', fpath, p1, p2, num_poses[l1], num_poses[l2]
                        os.system('rm {}'.format(fpath))
        except Exception as e:
            print 'mcss parse error', fpath
            print e
            print i,line
            os.system('rm {}'.format(fpath))

    return mcss_scores

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


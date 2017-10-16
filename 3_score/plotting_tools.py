import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as CM

def print_table(scores):
    #str1 = '|'
    #str2 = '|'
    #for i in ['min','ave','glide','opt']:# a.keys():#['min','ave','norm','opt','glide','us']:
    #    str1 = '{} {} |'.format(str1, i)
    #    str2 = '{} {} |'.format(str2, str(np.mean(a[i]))[:4])

    all_stats = {}

    success = lambda x: len([1 for i in x if i <= 2])

    top_str = '| '
    bot_str = '| '

    for n in [1, 5]:
        for f_name, f in {'ave':np.mean, 'med':np.median, 'suc':success}.items():
            for i in ['us', 'glide']:
                name = '{}-{}-{}'.format(i[0], n, f_name)
                val = scores.get_stats(i == 'glide', n, f)
                all_stats[name] = val
                top_str += name + ' | '
                bot_str += str(val)[:4] + ' | '
    print top_str
    print bot_str
    return all_stats
    ''' 
    glide_ave = str(np.mean(a['glide']))[:4]
    glide_med = str(np.median(a['glide']))[:4]
    us_med = str(np.median(a['opt']))[:4]
    us_ave = str(np.mean(a['opt']))[:4]
    us_success = str(sum([1.0/len(a['opt']) for i in a['opt'] if i <= 2.0]))[:4]
    glide_success = str(sum([1.0/len(a['glide']) for i in a['glide'] if i <= 2.0]))[:4]
    print '| min | ave | us (median) | glide (median) | us (ave) | glide (ave) | us (<2) | glide (<2)|'
    print '|{}|{}|{}|{}|{}|{}|{}|{}|'.format(str(np.mean(a['min']))[:4], str(np.mean(a['ave']))[:4],
                                             us_med, glide_med, us_ave, glide_ave, us_success, glide_success)
    #print str2'''

def plot_docking_by_structure(ligs, s, glides, n_list=[0,5,25,300], title='', scores=None):
    success_rate = {}
    plt.plot([2,2],[0,1],'--k')
    li = [l for l in ligs]# if s in glides[l]]
    prop_ligands = np.cumsum([-1.0/len(li)] + [1.0/len(li) for l in li] + [1.0/len(li)])
    for n in n_list:
        rmsds = [min([glides[l].poses[i].rmsd for i in range(min(n, max(glides[l].poses.keys()))+1)]) for l in li]
        success_rate[n] = sum([1 for i in rmsds if i <= 2])/float(len(rmsds))
        rmsds.sort()
        rmsds = [0] + rmsds + [rmsds[-1]]

        plt.step(rmsds, prop_ligands, label=n)
        
    if scores is not None: # len(us_rmsds) > 0:
        us_rmsds_1 = sorted([scores.get_top_rmsd(l, 1, False) for l in ligs])
        us_rmsds_5 = sorted([scores.get_top_rmsd(l, 5, False) for l in ligs])
        us_rmsds_1 = [0] + us_rmsds_1 + [us_rmsds_1[-1]]
        us_rmsds_5 = [0] + us_rmsds_5 + [us_rmsds_5[-1]]
        plt.step(us_rmsds_5, prop_ligands, label='us-5')
        plt.step(us_rmsds_1, prop_ligands, label='us-1')

    plt.gca().set_xlim([0,6])
    plt.gca().set_ylim([0,1])
    plt.xlabel('RMSD [A]')
    plt.ylabel('Cumulative Proportion of Ligands')
    plt.title(title+' '+s)
    plt.legend()
    plt.show()
    return success_rate
    
def plot_docking_output(ligs, struct, xglides, n, gscore=True):
    for lig in ligs:
        print lig

        if struct not in xglides[lig]: 
            print 'struct not found. struct, lig: ', struct, lig
            continue
        poses = xglides[lig][struct].poses
        pnum = [i for i in poses.keys()]
        g1 = [poses[i].gscore for i in poses.keys()]
        rmsd = [poses[i].rmsd for i in poses.keys()]
        if gscore:
            plt.plot(g1[:n], rmsd[:n], '.', markersize=10, label=lig)
        else:
            plt.plot(pnum[:n], rmsd[:n], '.', markersize=10, label=lig)

    #plt.gca().set_ylim([0,12])
    #plt.gca().set_xlim([0,25])
    plt.legend()
    plt.show()

def plot_all_poses(l1, l2, scores, lab, scores2=None, lab2=None):
    # plot pair score vs rmsd ave with arrows connecting pair 1 to pair 2
    # if a second set of scores is given, outputs a plot comparing the two scores
    
    # -1: do not plot the crystal pair scored against anything else
    
    all_scores = scores.get_all_scores(l1,l2)
    c_score = all_scores[-1][-1]
    plt.plot([0],[c_score], marker='*',markersize=10,color='red',label=lab)

    if scores2 is not None:
        all_scores2 = scores2.get_all_scores(l1,l2)
        c_score2 = all_scores2[-1][-1]
        plt.plot([0],[c_score2], marker='*',markersize=10,color='blue',label=lab2)
        plt.arrow(0,c_score,0,c_score2-c_score,fc='k',ec='k')

    rmsds1 = scores.get_rmsds(l1)
    rmsds2 = scores.get_rmsds(l2)
    
    points = []

    for p1 in range(np.shape(all_scores)[0] - 1):
        points = []
        for p2 in range(np.shape(all_scores)[1] - 1):
            x = (rmsds1[p1] + rmsds2[p2])/2.0
            y = all_scores[p1][p2]
            points.append((x,y))
            #plt.plot([x],[y],'r.')

            #if scores2 is not None:
            #    y2 = all_scores2[p1][p2]
            #    plt.plot([x],[y2],'b.')
            #    plt.arrow(x,y,0,y2-y,fc='k',ec='k')

        points.sort(key=lambda (x,y): y)
        points.reverse()
        #top = int(len(points)**0.5)
        x = [p[0] for p in points]
        y = [p[1] for p in points]
        plt.plot(x, y, '.')
        #plt.plot(x[top:], y[top:], 'b.')

    plt.legend()
    plt.xlabel('(rmsd1 + rmsd2)/2')
    plt.ylabel('pair score')
    plt.gca().set_xlim([0,12])
    #plt.gca().set_ylim([0,700])
    plt.show()

def plot_clusters(scores, lab=''):
    for l in scores.ligands:
        for p in scores.pose_clusters[l]:
            if p == -1: symb = 'k*'
            else: symb = 'k.'
            rmsd = np.mean([scores.get_rmsds(li)[pi] for (li,pi) in scores.pose_clusters[l][p]])
            plt.plot([rmsd],[scores.objective(scores.pose_clusters[l][p])], symb)
    plt.show()

def plot_magnitudes(scores, title=''):
    fig, ax = plt.subplots()

    top_poses = [np.argmax(scores.get_final_scores(l)[:-1]) for l in scores.ligands]
    top_pose_magnitudes = [0 for l in scores.ligands]    

    highest_magnitude_poses = {}
    for (i, l) in enumerate(scores.ligands):

        fp_top = scores.glides[l][scores.struct].poses[top_poses[i]].fp
        top_pose_magnitudes[i] = scores.score_pose_pair(fp_top, fp_top)

        highest_magnitude_poses[l] = (-1, 0)
        for p in range(scores.num_poses[l]):
            fp = scores.glides[l][scores.struct].poses[p].fp
            norm = scores.score_pose_pair(fp, fp)
            if norm > highest_magnitude_poses[l][1]:
                highest_magnitude_poses[l] = (p, norm)
    
    max_pose_magnitudes = [highest_magnitude_poses[l][1] for l in scores.ligands]
    plt.plot(max_pose_magnitudes, marker='.', markersize=10, color='red', label='max norm')
    plt.plot(top_pose_magnitudes, marker='.', markersize=10, color='blue',label='top norm')
    
    plt.legend()
    ax.set_xticklabels(scores.ligands, minor=False, rotation='vertical')
    ax.set_xticks(np.arange(0,len(scores.ligands),1))

    ax.set_title(title)
    ax.set_xlabel('Ligands')
    ax.set_ylabel('Magnitude')

    plt.show()

    max_norm_rmsds = [scores.get_rmsds(l)[highest_magnitude_poses[l][0]] for l in scores.ligands]
    plot_final_rmsds(scores, title=title, max_norm_rmsds=max_norm_rmsds)

def get_interactions(pose_cluster, interactions):
    all_residues = []
    for (l,p), fp in pose_cluster.items():
        all_residues.extend([r for r in fp if r not in all_residues])
    
    shared = {}
    empty_fp = [0 for i in range(max(interactions)+1)]
    for r in all_residues:
        for i in interactions:
            strengths = []
            for j in range(len(pose_cluster.keys())):
                for k in range(j+1, len(pose_cluster.keys())):
                    fp1 = pose_cluster[pose_cluster.keys()[j]]
                    fp2 = pose_cluster[pose_cluster.keys()[k]]
                    strengths.append(fp1.get(r,empty_fp)[i]*fp2.get(r,empty_fp)[i])
            if max(strengths) > 0:
                shared[(r,i)] = sum(np.sqrt(strengths))/float(len(pose_cluster)*(len(pose_cluster) - 1))

    return shared, sum([s for k, s in shared.items()]) # sorted(shared.keys(), key=lambda x:-shared[x])[:max_r]

def plot_shared_interactions(scores, max_r=10, title=''):
    interactions = [i for i in range(12)]
    
    clusters = {'crystal':{(l,-1):scores.crystals[l] for l in scores.ligands}}
    for c in ['opt', 'glide']:
        clusters[c] = {}
        for i, l in enumerate(scores.ligands):
            pnum = scores.all_analysis[c][0][i]
            clusters[c][(l, pnum)] = scores.glides[l].poses[pnum].fp
    
    # pose cluster is a dict mapping (lig_name, pnum) to fp instance
    shared, s1 = get_interactions(clusters['crystal'], interactions)
    to_plot = sorted(shared.keys(), key=lambda x:-shared[x])[:max_r]

    shared2, s2 = get_interactions(clusters['opt'], interactions)
    to_plot2 = sorted(shared2.keys(), key=lambda x:-shared2[x])[:max_r]
    
    shared3, s3 = get_interactions(clusters['glide'], interactions)
    to_plot3 = sorted(shared3.keys(), key=lambda x:-shared3[x])[:max_r]
    
    to_plot.extend([i for i in to_plot2 if i not in to_plot])
    to_plot.extend([i for i in to_plot3 if i not in to_plot])
    
    to_plot.sort(key=lambda x: -(shared.get(x,0) + shared2.get(x,0) + shared3.get(x,0)))
    
    plt.plot([shared3.get(i, 0) for i in to_plot], '.-', markersize=10, label='glide score: {}'.format(s3))
    plt.plot([shared2.get(i, 0) for i in to_plot], '.-', markersize=10, label='opt score: {}'.format(s2))
    plt.plot([shared.get(i, 0) for i in to_plot], '.-', markersize=10, label='crystal score: {}'.format(s1))

    plt.gca().set_xticklabels(to_plot, minor=False, rotation='vertical')
    plt.gca().set_xticks(np.arange(0, len(to_plot), 1))
    plt.title(title)
    plt.ylim([-0.05, 1.05])
    plt.legend()
    plt.show()

def plot_final_rmsds(scores, title='', scores2=None, lab2='', max_norm_rmsds=None):#, show_glide=False):
    
    a = scores.all_analysis
    for i in ['min', 'ave', 'glide', 'opt']:# 'norm', 'opt', 'glide', 'us']:
        plt.plot(a[i][1][:], marker='.', markersize=10, label='{}: {}'.format(i, str(np.mean(a[i][1][:]))[:4]))

    if scores2 is not None:
        if False in [scores.ligands[i] == scores2.ligands[i] for i in range(len(scores.ligands))]:
            raise Exception

        us_rmsds2 = scores2.all_analysis['us'][1][:]
        plt.plot(us_rmsds2, marker='.', markersize=10, color='blue',label=lab2+str(np.mean(us_rmsds2))[:4])

    plt.legend()
    plt.gca().set_xticklabels(scores.ligands, minor=False, rotation='vertical')
    plt.gca().set_xticks(np.arange(0,len(scores.ligands),1))
    plt.gca().set_ylim([0,12])
    
    plt.gca().set_title(title)
    plt.gca().set_xlabel('Ligands')
    plt.gca().set_ylabel('RMSD [A]')

    plt.show()

def plot_final_scores(scores, lab=''):
    top_pose_scores = [np.max(scores.get_final_scores(l)[:-1]) for l in scores.ligands]
    crystal_pose_scores = [scores.get_final_scores(l)[-1] for l in scores.ligands]

    fig, ax = plt.subplots()
    plt.plot(crystal_pose_scores, marker='*', markersize=10, color="red")
    plt.plot(top_pose_scores, marker='.', markersize=10, color="red")
    ax.set_xticks(np.arange(0,len(scores.ligands),1))
    ax.set_xticklabels(scores.ligands, minor=False, rotation='vertical')
    plt.show()

def plot_scores_vs_rmsds(l, scores, lab='', scores2=None, lab2=''):
    final_scores = scores.get_final_scores(l)
    rmsds = scores.get_rmsds(l)
    fig, ax = plt.subplots()
    plt.plot(rmsds[:-1], final_scores[:-1], marker='.', markersize=10, color='red', label=lab, linestyle='None')
    plt.plot([0], [final_scores[-1]], marker='*', markersize=10, color='red')

    if scores2 is not None:
        final_scores2 = scores2.get_final_scores(l)
        rmsds2 = scores2.get_rmsds(l)
        plt.plot(rmsds2[:-1], final_scores2[:-1], marker='.', markersize=10, color='blue', label=lab2, linestyle='None')
        plt.plot([0], [final_scores2[-1]], marker='*', markersize=10, color='blue')
    plt.legend()
    plt.gca().set_xlim([0,12])
    plt.show()

def plot_ifp_comparison(fp1, fp2, i=[0,1,2,3], lab1='fp1', lab2='fp2', title={0:'res=hdonor',1:'res=hacc',2:'sb',3:'lj'}):

    is_interaction = lambda fp, r, i: r in fp and True in [fp[r][j] >= 1 for j in i]

    sort_key = lambda r: int(r)
    fp1_res = sorted([r for r in fp1 if is_interaction(fp1, r, i) and not is_interaction(fp2, r, i)], key=sort_key)
    both_res= sorted([r for r in fp1 if is_interaction(fp1, r, i) and is_interaction(fp2, r, i)], key=sort_key)
    fp2_res = sorted([r for r in fp2 if not is_interaction(fp1, r, i) and is_interaction(fp2, r, i)], key=sort_key)

    symb = {0:'*',1:'s',2:'d',3:'.'}

    res = fp1_res + both_res + fp2_res
    val1 = {j:[] for j in i}
    val2 = {j:[] for j in i}

    for j in i:
        for r in res:
            if is_interaction(fp1, r, [j]):
                val1[j].append(fp1[r][j])
            else:
                val1[j].append(-1)
            if is_interaction(fp2, r, [j]):
                val2[j].append(fp2[r][j])
            else:
                val2[j].append(-1)

        plt.plot(val1[j], symb[j]+'r', markersize=10, label=title[j])
        plt.plot(val2[j], symb[j]+'b', markersize=10)

    plt.legend()
    plt.title('red = {}, blue = {}'.format(lab1, lab2))
    plt.gca().set_ylim([0,12])
    plt.gca().set_xticklabels(res, minor=False, rotation='vertical')
    plt.gca().set_xticks(np.arange(0,len(res),1))
    plt.show()

def plot_score_breakdown(scores, lab=''):
    final_poses = [np.argmax(scores.get_final_scores(l)[:-1]) for l in scores.ligands]
    bd = {l: scores.score_breakdown(l, final_poses[i]) for i, l in enumerate(scores.ligands)}

    lig_ind = np.arange(len(scores.ligands))

    leg = {0:'res=hdonor', 1:'res=hacc', 2:'sb', 3:'lj', 4:''}

    plt.bar(lig_ind, [bd[l][0] for l in scores.ligands], label=leg[0])
    for i in range(1,len(bd[scores.ligands[0]])):
        prev = [sum([bd[l][i-n-1] for n in range(i)]) for l in scores.ligands]
        plt.bar(lig_ind, [bd[l][i] for l in scores.ligands], bottom=prev, label=leg[i])
    plt.gca().set_xticklabels(scores.ligands, minor=False, rotation='vertical')
    plt.gca().set_xticks(np.arange(0,len(scores.ligands),1))
    plt.legend()
    plt.show()

def plot_n_rmsds(scores, lab=''):
    final_poses = [np.argmax(scores.get_final_scores(l)[:-1]) for l in scores.ligands]
    n_rmsds = [scores.rmsd_of_neighbors(l, final_poses[i]) for i, l in enumerate(scores.ligands)]
    plt.plot(n_rmsds, marker='.', markersize=10)
    plt.gca().set_xticklabels(scores.ligands, minor=False, rotation='vertical')
    plt.gca().set_xticks(np.arange(0,len(scores.ligands),1))
    plt.gca().set_ylim([0,12])
    plt.show()


def heatmap(A, row_vals, col_vals, red=10.25, row_label='ligands', col_label='structures'):
    fig, ax = plt.subplots()

    cmap = CM.jet
    cmap.set_under('black')

    heatmap = ax.pcolor(A, cmap=cmap, vmin=0, vmax=red)
    plt.colorbar(heatmap)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(A.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(A.shape[0]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    #labels
    #column_labels = x_vals#ligstructs
    #row_labels = y_vals#gridstructs
    ax.set_xticklabels(col_vals, minor=False, rotation = 'vertical')
    ax.set_yticklabels(row_vals, minor=False)
    ax.plot(ax.get_xlim(), ax.get_ylim()[::-1], linewidth = 4, c="m")
    plt.title('')
    plt.xlabel(col_label)
    plt.ylabel(row_label)
    #plt.savefig(plotTitle + '.png')
    plt.show()

def refine_poses(glides, require_fp=False):
    # glides maps [lig][struct] to ligands
    # goal: output poses. map [lig] to all poses, ordered by gscore
    poses = {}
    for lig in glides.keys():
        poses[lig] = []
        for struct in glides[lig].keys():
            if lig == struct: continue
            p = glides[lig][struct].poses
            poses[lig].extend([p[i] for i in p.keys() if (not require_fp) or (p[i].fp is not None)])
        poses[lig].sort(key=lambda x: x.gscore)
    return poses

def helpfully_frozen(ligs, structs, glides, n):
    luck = 0
    total = 0
    for l in ligs:
        for s in structs:
            if not (l in glides and s in glides[l]): continue
            total += 1
            p = [glides[l][s].poses[i].rmsd for i in range(n) if i < len(glides[l][s].poses)]
            if np.mean(p) < 2 and np.var(p) < 1:
                luck += 1
    print 'got lucky on {} of {} pairs.'.format(luck, total)

def get_docking_stats(ligs, structs, glides, n, func):
    A = np.zeros( (len(structs), len(ligs)) )
    for i, struct in enumerate(structs):
        for j, lig in enumerate(ligs):
            if lig in glides.keys() and struct in glides[lig].keys():
                num_poses = min(n, len(glides[lig][struct].poses.keys()))
                A[i, j] = func([glides[lig][struct].poses[k].rmsd for k in range(num_poses)])
            else: A[i, j] = np.nan
    return A

def show_results_for_all_ligand_pairs(scores, lab='', scores2=None, lab2=''):    
    count = 0
    pair_indices = {}
    ave_change = 0
    for i1 in range(len(scores.ligands)):
        for i2 in range(i1 + 1, len(scores.ligands)):
                        
            rmsds1 = scores.get_rmsds(scores.ligands[i1])[:-1]
            rmsds2 = scores.get_rmsds(scores.ligands[i2])[:-1]
            all_scores = scores.get_all_scores(scores.ligands[i1], scores.ligands[i2])
            (p1, p2) = np.shape(all_scores)
            
            # 1: lowest/highest rmsd pair available
            min_rmsd = (np.min(rmsds1) + np.min(rmsds2))/2.0
            max_rmsd = (np.max(rmsds1) + np.max(rmsds2))/2.0

            plt.plot([count], [min_rmsd], marker='d', markersize=25, color="black")
            plt.plot([count], [max_rmsd], marker='d', markersize=25, color="black")
            plt.arrow(count,min_rmsd,0,max_rmsd - min_rmsd, fc='k', ec='k')
            
            # 2: our performance 1
            pair = np.unravel_index(np.argmax(all_scores[:-1][:-1]), (p1 - 1, p2 - 1))
            top_score = all_scores[pair]
            top_rmsd = (rmsds1[pair[0]] + rmsds2[pair[1]])/2.0
            plt.plot([count], [top_rmsd], marker='.', markersize=25, color="red")
            
            pair_indices[count] = (scores.ligands[i1], scores.ligands[i2], top_rmsd, min_rmsd, pair)
            
            # 3: (optional) our performance 2
            if scores2 is not None:
                all_scores2 = scores2.get_all_scores(scores.ligands[i1], scores.ligands[i2])
            
                pair2 = np.unravel_index(np.argmax(all_scores2[:-1][:-1]), (p1-1, p2-1))
                top_rmsd2 = (rmsds1[pair2[0]] + rmsds2[pair2[1]])/2.0
                plt.plot([count], [top_rmsd2], marker='.', markersize=25, color="blue")
                ave_change += top_rmsd2 - top_rmsd
            
            count += 1    
    
    ave_change /= float(count)

    if scores2 is not None:
        print 'Average rmsd change = ' + str(ave_change)
        print 'If that number is negative, scores2 is better than scores1!'
        
    plt.arrow(-1,1,2+count,0, fc='k', ec='k')
    plt.arrow(-1,2,2+count,0, fc='k', ec='k')
    plt.arrow(-1,4,2+count,0, fc='k', ec='k')

    
    plt.ylabel('(rmsd1 + rmsd2)/2')
    plt.xlabel('ligand pair')
    plt.tick_params(axis='both', which='major', labelsize=30)
    
    #plt.legend()
    plt.show()
    return pair_indices

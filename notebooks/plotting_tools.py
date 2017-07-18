import numpy as np
from matplotlib import pyplot as plt
import matplotlib.cm as CM

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

    for p1 in range(np.shape(all_scores)[0] - 1):
        for p2 in range(np.shape(all_scores)[1] - 1):
            x = (rmsds1[p1] + rmsds2[p2])/2.0
            y = all_scores[p1][p2]
            plt.plot([x],[y],'r.')

            if scores2 is not None:
                y2 = all_scores2[p1][p2]
                plt.plot([x],[y2],'b.')
                plt.arrow(x,y,0,y2-y,fc='k',ec='k')

    plt.legend()
    plt.xlabel('(rmsd1 + rmsd2)/2')
    plt.ylabel('pair score')
    #plt.gca().set_xlim([-0.5,9])
    #plt.gca().set_ylim([0,700])
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
        

def plot_final_rmsds(scores, title='', scores2=None, lab2='', max_norm_rmsds=None):#, show_glide=False):
    fig, ax = plt.subplots()

    top_pose_rmsds = [scores.get_rmsds(l)[np.argmax(scores.get_final_scores(l)[:-1])] for l in scores.ligands]
    min_rmsds = [np.min(scores.get_rmsds(l)[:-1]) for l in scores.ligands]
    plt.plot(min_rmsds, marker='d', markersize=10, color="black", label='best: '+str(np.mean(min_rmsds))[:4])

    if max_norm_rmsds is not None:
        plt.plot(max_norm_rmsds, marker='.', markersize=10, color='magenta', label='max norm: '+str(np.mean(max_norm_rmsds))[:4])

    if True:#show_glide:
        glide_rmsds = [scores.get_rmsds(l)[np.argmin(scores.get_gscores(l))] for l in scores.ligands]
        plt.plot(glide_rmsds, marker='.', markersize=10, color='green', label='glide: '+str(np.mean(glide_rmsds))[:4])

    plt.plot(top_pose_rmsds, marker='.', markersize=10, color="red",label='us: '+str(np.mean(top_pose_rmsds))[:4])

    if scores2 is not None:
        top_rmsds2 = [scores2.get_rmsds(l)[np.argmax(scores2.get_final_scores(l)[:-1])] for l in scores2.ligands]
        plt.plot(top_rmsds2, marker='.', markersize=10, color='blue',label=lab2+str(np.mean(top_rmsds2))[:4])

    plt.legend()
    ax.set_xticklabels(scores.ligands, minor=False, rotation='vertical')
    ax.set_xticks(np.arange(0,len(scores.ligands),1))
    plt.gca().set_ylim([0,10])
    
    ax.set_title(title)
    ax.set_xlabel('Ligands')
    ax.set_ylabel('RMSD [A]')

    plt.show()

    return [(scores.ligands[i],top_pose_rmsds[i],min_rmsds[i],glide_rmsds[i]) for i in range(len(scores.ligands))]

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
    plt.plot(rmsds[:-1], final_scores[:-1], marker='.', markersize=10, color='red', label=lab, linestyle='None')
    plt.plot([0], [final_scores[-1]], marker='*', markersize=10, color='red')

    if scores2 is not None:
        final_scores2 = scores2.get_final_scores(l)
        plt.plot(rmsds[:-1], final_scores2[:-1], marker='.', markersize=10, color='blue', label=lab2, linestyle='None')
        plt.plot([0], [final_scores2[-1]], marker='*', markersize=10, color='blue')
    plt.legend()
    plt.show()

def heatmap(A, glides):
    fig, ax = plt.subplots()

    cmap = CM.jet
    cmap.set_under('black')

    heatmap = ax.pcolor(A, cmap=cmap, vmin=0, vmax=10.25)
    plt.colorbar(heatmap)

    # put the major ticks at the middle of each cell
    ax.set_xticks(np.arange(A.shape[1]) + 0.5, minor=False)
    ax.set_yticks(np.arange(A.shape[0]) + 0.5, minor=False)
    ax.invert_yaxis()
    ax.xaxis.tick_top()

    #labels
    column_labels = glides.keys()
    row_labels = glides.keys()
    ax.set_xticklabels(column_labels, minor=False, rotation = 'vertical')
    ax.set_yticklabels(row_labels, minor=False)
    ax.plot(ax.get_xlim(), ax.get_ylim()[::-1], linewidth = 4, c="m")
    plt.title('')
    #plt.savefig(plotTitle + '.png')
    plt.show()

 # Top scoring GLIDE pose
def top_pose(glides):
    A = np.zeros( (len(glides.keys()), len(glides.keys())) )
    for i, grid in enumerate(glides):
        for j, ligand in enumerate(glides):
            A[i, j] = glides[ligand][grid].poses[0].rmsd
    return A

 # Top scoring GLIDE pose in top n poses
def best_pose(ligstructs, gridstructs, glides, n):
    A = np.zeros( (len(glides.keys()), len(glides.keys())) )
    for i, grid in enumerate(gridstructs):
        for j, lig in enumerate(ligstructs):
            if lig in glides.keys() and grid in glides[lig].keys():
                num_poses = min(n, len(glides[lig][grid].poses.keys()))
                #if num_poses == 0: A[i, j] = np.nan
                A[i, j] = min([glides[lig][grid].poses[k].rmsd for k in range(num_poses)])
            else: A[i, j] = np.nan
    return A

def get_structure_and_ligands(rmsd_matrix, rmsd_filter, structures):
    (best_struct, lig, score) = (None, [], 10)
    for i, struct in enumerate(structures):
        good_lig = [(l, rmsd_matrix[i][j]) for j, l in enumerate(structures) if rmsd_filter(rmsd_matrix[i][j])]
        ave_rmsd = sum([rmsd for (l, rmsd) in good_lig])/float(len(good_lig))
        if len(good_lig) > len(lig) or (len(good_lig) == len(lig) and ave_rmsd < score):
            (best_struct, lig, score) = (struct, good_lig, ave_rmsd)
    
    return (best_struct, [i[0] for i in lig], score)

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

import numpy as np

def objective_partial(pose_cluster, k1, phys_const, pair_const):
        score = phys_const*pose_cluster[k1].gscore
        for k, p in pose_cluster.items():
            if k != k1: 
                score += pair_const*score_fp_pair(pose_cluster[k1].fp, p.fp)/float(len(pose_cluster.keys()))
        return score

def objective(pose_cluster, phys_const, pair_const):
    score = 0
    ligs = pose_cluster.keys()
    n = float(len(ligs))
    for i, l1 in enumerate(ligs):
        score += phys_const*pose_cluster[ligs[i]].gscore/n
        for j in range(i + 1, len(ligs)):
            l2 = ligs[j]
            score += pair_const*score_fp_pair(pose_cluster[ligs[i]].fp, pose_cluster[ligs[j]].fp)/n**2
    return score

def score_fp_pair(fp1, fp2):
    score = 0
    for r in fp1:
        if r in fp2:
            score += np.dot(np.sqrt(fp1[r]), np.sqrt(fp2[r]))
    return score


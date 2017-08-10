import numpy as np
from fingerprint import FuzzyFingerPrint

class Scores:
    def __init__(self, glides, crystals, ligands, structure, n):
        self.all_rmsds = {}
        self.all_gscores = {}

        self.glides = glides
        self.crystals = crystals
        self.ligands = [l for l in ligands if structure in glides[l]]
        self.struct = structure
        self.num_poses = {l: min(n,len(glides[l][structure].poses)) for l in self.ligands}
        # the crystal poses are the last row/column of each matrix here
        self.all_scores = self.score_all_pairs_of_poses()

        self.optimized_scores = self.joint_optimize()
        self.all_analysis = self.aggregate_scores()

    def aggregate_scores(self):
        all_analysis = {
            'opt'   : [self.optimized_scores[l] for l in self.ligands],
            'glide' : [0 for l in self.ligands],
            'min'   : [np.argmin(self.get_rmsds(l)[:-1]) for l in self.ligands],
        }
        
        for i in all_analysis:
            t = np.zeros( (2, len(self.ligands)) )
            for j, p in enumerate(all_analysis[i]):
                t[0][j] = p
                t[1][j] = self.get_rmsds(self.ligands[j])[p]
            all_analysis[i] = t
        all_analysis['ave'] = np.zeros( (2, len(self.ligands)) )
        all_analysis['ave'][0][:] = [None for i in self.ligands]
        all_analysis['ave'][1][:] = [np.mean(self.get_rmsds(l)[:-1]) for l in self.ligands]
        
        return all_analysis

    def objective_partial(self, pose_cluster, l1):
        score = 0
        for l2 in self.ligands:
            if l1 == l2: continue
            score += self.score_pose_pair(l1, pose_cluster[l1], l2, pose_cluster[l2])
        return score

    def objective(self, pose_cluster):
        score = 0
        for i in range(len(self.ligands)):
            for j in range(i + 1, len(self.ligands)):
                l1, l2 = self.ligands[i], self.ligands[j]
                score += self.score_pose_pair(l1, pose_cluster[l1], l2, pose_cluster[l2])
        return score

    def joint_optimize(self):
        score, cluster = self.optimize_helper({l: 0 for l in self.ligands})
        for i in range(5):
            new_score, new_cluster = self.optimize_helper({l: np.random.randint(self.num_poses[l]) for l in self.ligands})
            if new_score > score:
                score, cluster = new_score, new_cluster
        return cluster

    def optimize_helper(self, init_cluster):
        # randomly pick a ligand
        # optimize that ligand
        time_since_update = 0
        pose_cluster = init_cluster # {l: 0 for l in self.ligands}
        while time_since_update < len(self.ligands)**2:
            time_since_update += 1

            rand_lig = np.random.choice(self.ligands)
            old_rmsd = self.get_rmsds(rand_lig)[pose_cluster[rand_lig]]
    
            (max_score, pose_num) = (0, 0)
            for p in range(self.num_poses[rand_lig]):
                pose_cluster[rand_lig] = p
                new_score = self.objective_partial(pose_cluster, rand_lig)
                if new_score > max_score:
                    max_score, pose_num = new_score, p
            pose_cluster[rand_lig] = pose_num
    
            new_rmsd = self.get_rmsds(rand_lig)[pose_cluster[rand_lig]]
    
            if new_rmsd != old_rmsd:
                time_since_update = 0
        return self.objective(pose_cluster), pose_cluster

    def score_pose_pair(self, l1, p1, l2, p2):
        if hasattr(self, 'all_scores') and (l1, l2) in self.all_scores:
            return self.all_scores[(l1, l2)][p1][p2]
        if hasattr(self, 'all_scores') and (l2, l1) in self.all_scores:
            return self.all_scores[(l2, l1)][p2][p1]

        fp1 = self.glides[l1][self.struct].poses.get(p1, self.crystals[l1]).fp
        fp2 = self.glides[l2][self.struct].poses.get(p2, self.crystals[l2]).fp
        score = 0
        for r in fp1.feats:
            if r in fp2.feats:
                score += np.dot(fp1.feats[r], fp2.feats[r])
        return score

    def score_all_pairs_of_poses(self):
        all_scores = {}
        for i1 in range(len(self.ligands)):
            for i2 in range(i1+1, len(self.ligands)):            
                l1, l2 = self.ligands[i1], self.ligands[i2]
                d1, d2 = self.num_poses[l1], self.num_poses[l2]
                # +-1 allows an extra (last) row/column to be used to score the crystallized poses
                all_scores[(l1, l2)] = np.zeros((d1 + 1, d2 + 1)) 
                
                for p1 in range(-1,d1):
                    for p2 in range(-1,d2):
                        all_scores[(l1, l2)][p1][p2] = self.score_pose_pair(l1, p1, l2, p2)
        return all_scores

    def get_best_cluster(self):
        cluster = {}
        for l in self.ligands:
            best_p = int(np.argmin(self.get_rmsds(l)[:-1]))
            cluster[(l,best_p)] = self.glides[l][self.struct].poses[best_p].fp
        return cluster

    def get_rmsds(self,l):
        # returns all rmsds and then 0 for the crystal pose
        if l not in self.all_rmsds: 
            self.all_rmsds[l] = [self.glides[l][self.struct].poses[p].rmsd for p in range(self.num_poses[l])] + [0]
        return self.all_rmsds[l]

    def get_gscores(self,l):
        if l not in self.all_gscores:
            self.all_gscores[l] = [self.glides[l][self.struct].poses[p].gscore for p in range(self.num_poses[l])]
        return self.all_gscores[l]

import numpy as np

from objective import objective, objective_partial, score_fp_pair

class Scores:
    def __init__(self, glides, crystals, ligands, structure, n, gscore_weight=0, all_scores=None):
        self.all_rmsds = {}
        self.all_gscores = {}

        self.gscore_weight = -1*gscore_weight

        self.glides = glides
        self.crystals = crystals
        self.ligands = ligands#[l for l in ligands if structure in glides[l]]
        self.struct = structure
        
        self.num_poses = {l: min(n,len(glides[l].poses)) for l in self.ligands}
        # the crystal poses are the last row/column of each matrix here
        if all_scores:
            self.all_scores = all_scores
        else:
            self.all_scores = self.score_all_pairs_of_poses()

        self.optimized_scores = self.joint_optimize()
        self.all_analysis = self.aggregate_scores()

    def aggregate_scores(self):
        all_analysis = {
            'opt'   : [self.optimized_scores[l] for l in self.ligands],
            'glide' : [0 for l in self.ligands],
            'min'   : [np.argmin(self.get_rmsds(l)) for l in self.ligands],
        }
        
        for i in all_analysis:
            t = np.zeros( (2, len(self.ligands)) )
            for j, p in enumerate(all_analysis[i]):
                t[0][j] = p
                t[1][j] = self.get_rmsds(self.ligands[j])[p]
            all_analysis[i] = t
        all_analysis['ave'] = np.zeros( (2, len(self.ligands)) )
        all_analysis['ave'][0][:] = [None for i in self.ligands]
        all_analysis['ave'][1][:] = [np.mean(self.get_rmsds(l)) for l in self.ligands]
        
        return all_analysis

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
        pose_cluster = init_cluster

        phys_const = self.gscore_weight
        pair_const = 1

        while time_since_update < len(self.ligands)**2:
            time_since_update += 1

            rand_lig = np.random.choice(self.ligands)
            
            (max_score, pose_num) = (objective_partial(pose_cluster, rand_lig, phys_const, pair_const), pose_cluster[rand_lig])
            for p in range(self.num_poses[rand_lig]):
                pose_cluster[rand_lig] = p
                new_score = objective_partial(pose_cluster, rand_lig, phys_const, pair_const)
                if new_score > max_score:
                    max_score, pose_num = new_score, p
                    time_since_update = 0
            pose_cluster[rand_lig] = pose_num
    
        return objective(pose_cluster, phys_const, pair_const), pose_cluster

    def score_pose_pair(self, l1, p1, l2, p2):
        if hasattr(self, 'all_scores') and (l1, l2) in self.all_scores:
            return self.all_scores[(l1, l2)][p1][p2]
        if hasattr(self, 'all_scores') and (l2, l1) in self.all_scores:
            return self.all_scores[(l2, l1)][p2][p1]

        if p1 == -1:
            fp1 = self.crystals[l1]
        else:
            fp1 = self.glides[l1].poses[p1].fp
        if p2 == -1:
            fp2 = self.crystals[l2]
        else:
            fp2 = self.glides[l2].poses[p2].fp

        return score_fp_pair(fp1, fp2)

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
            best_p = int(np.argmin(self.get_rmsds(l)))
            cluster[(l,best_p)] = self.glides[l].poses[best_p].fp
        return cluster

    def get_rmsds(self,l):
        # returns all rmsds and then 0 for the crystal pose
        if l not in self.all_rmsds: 
            self.all_rmsds[l] = [self.glides[l].poses[p].rmsd for p in range(self.num_poses[l])] # + [0]
        return self.all_rmsds[l]

    def get_gscores(self,l):
        if l not in self.all_gscores:
            self.all_gscores[l] = [self.glides[l].poses[p].gscore for p in range(self.num_poses[l])]
        return self.all_gscores[l]

    def get_scores(self, l):
        if not hasattr(self, 'all_opt_scores') or l not in self.all_opt_scores:
            self.all_opt_scores = {}
            self.all_opt_scores[l] = []
            for p in range(self.num_poses[l]):
                new_cluster = {i:self.optimized_scores[i] for i in self.optimized_scores}
                new_cluster[l] = p
                self.all_opt_scores[l].append(objective(new_cluster, self.gscore_weight, 1))
        return self.all_opt_scores[l]

    def get_stats(self, glide, n, func):
        return func([self.get_top_rmsd(l, n, glide) for l in self.ligands])
        '''
        all_stats = {}
        all_stats[('us', 1, 'median')] = np.median([self.get_top(l, 1, False) for l in self.ligands])
        all_stats[('us', 5, 'median')] = np.median([self.get_top(l, 5, False) for l in self.ligands])
        all_stats[('us', 1, 'ave')] = np.mean([self.get_top(l, 1, False) for l in self.ligands])
        all_stats[('us', 5, 'ave')] = np.mean([self.get_top(l, 5, False) for l in self.ligands])
        all_stats[('us', 1-success')] = len([1 for l in self.ligands if self.get_top(l, 1, False) <= 2])
        all_stats[('us-5-success')] = len([1 for l in self.ligands if self.get_top(l, 5, False) <= 2])
        all_stats[('gl-1-median')] = np.median([self.get_top(l, 1, True ) for l in self.ligands])
        all_stats[('gl-5-median')] = np.median([self.get_top(l, 5, True ) for l in self.ligands])
        all_stats[('gl-1-ave')] = np.mean([self.get_top(l, 1, True) for l in self.ligands])
        all_stats[('gl-5-ave')] = np.mean([self.get_top(l, 5, True) for l in self.ligands])
        all_stats[('gl-1-success')] = len([1 for l in self.ligands if self.get_top(l, 1, True) <= 2])
        all_stats[('gl-5-success')] = len([1 for l in self.ligands if self.get_top(l, 5, True) <= 2])
        return all_stats'''

    def get_top_num(self, l, n, glide=True):
        if glide:
            return np.argmin([self.glides[l].poses[p].rmsd for p in range(n)])
        our_scores = [(i, self.get_scores(l)[i]) for i in range(self.num_poses[l])]
        our_scores.sort(key=lambda x: -x[1])
        best_pose = np.argmin([self.glides[l].poses[our_scores[i][0]].rmsd for i in range(n)])
        return our_scores[best_pose][0]

    def get_top_rmsd(self, l, n, glide=True):
        if glide:
            return np.min([self.glides[l].poses[p].rmsd for p in range(n)])
        our_scores = [(i, self.get_scores(l)[i]) for i in range(self.num_poses[l])]
        our_scores.sort(key=lambda x: -x[1])
        return np.min([self.glides[l].poses[our_scores[i][0]].rmsd for i in range(n)])

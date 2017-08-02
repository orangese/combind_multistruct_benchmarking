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
        self.scores_from_one_ligand = self.rank_poses_based_on_one_other_ligand()
        self.final_scores_for_each_ligand = self.rank_poses_based_on_all_other_ligands()

        self.optimized_scores = self.joint_optimize()

        self.all_analysis = self.aggregate_scores()

    def norm(self, l, p):
        fp = self.glides[l][self.struct].poses[p].fp
        return self.score_pose_pair(fp, fp)

    def aggregate_scores(self):
        norms = {l:[self.norm(l,p) for p in range(self.num_poses[l])] for l in self.ligands}

        all_analysis = {
            'us'    : [np.argmax(self.final_scores_for_each_ligand[l][:-1]) for l in self.ligands],
            'opt'   : [self.optimized_scores[l] for l in self.ligands],
            'glide' : [0 for l in self.ligands],
            'min'   : [np.argmin(self.get_rmsds(l)[:-1]) for l in self.ligands],
            'norm'  : [np.argmax(norms[l]) for l in self.ligands]
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

    def add_overlap(self, fp1, fp2, i=None):
        new_feats = {}
        for r in fp1.feats:
            new_feats[r] = np.add(fp1.feats[r], fp2.feats.get(r,[0,0,0,0,0]))
        for r in fp2.feats:
            if r not in fp1.feats:
                new_feats[r] = fp2.feats[r]
        return FuzzyFingerPrint(new_feats)

 
    def add_objective(self, pose_cluster):
        score = 0
        (l1, p1) = pose_cluster[0]
        (l2, p2) = pose_cluster[1]
        fp = self.add_overlap(self.glides[l1][self.struct].poses[p1].fp, self.glides[l2][self.struct].poses[p2].fp)
        for li, pi in pose_cluster[2:]:
            fp = self.add_overlap(fp, self.glides[l1][self.struct].poses[pi].fp)
        
        return np.sum([np.sum(np.square(fp.feats[r])) for r in fp.feats])


    def objective(self, pose_cluster):
        score = 0
        for l1, p1 in pose_cluster:
            for l2, p2 in pose_cluster:
                if l1 == l2: continue
                if p1 == -1: fp1 = self.crystals[l1].fp
                else: fp1 = self.glides[l1][self.struct].poses[p1].fp
                if p2 == -1: fp2 = self.crystals[l2].fp
                else: fp2 = self.glides[l2][self.struct].poses[p2].fp
                score += self.score_pose_pair(fp1, fp2)
        return score

    def joint_optimize(self):
        # randomly pick a ligand
        # optimize that ligand
        time_since_update = 0
        pose_cluster = [(l, np.argmax(self.get_final_scores(l)[:-1])) for l in self.ligands]
        while time_since_update < len(self.ligands)**2:
            time_since_update += 1

            rand_lig = np.random.choice(self.ligands)
            lig_ind = self.ligands.index(rand_lig)
    
            old_rmsd = self.get_rmsds(rand_lig)[pose_cluster[lig_ind][1]]
    
            (max_score, pose_num) = (0, 0)
            for p in range(self.num_poses[rand_lig]):
                pose_cluster[lig_ind] = (rand_lig, p)
                new_score = self.add_objective(pose_cluster)
                if new_score > max_score:
                    max_score, pose_num = new_score, p
            pose_cluster[lig_ind] = (rand_lig, pose_num)
    
            new_rmsd = self.get_rmsds(rand_lig)[pose_cluster[lig_ind][1]]
    
            if new_rmsd != old_rmsd:
                time_since_update = 0
        return {l: p for (l, p) in pose_cluster}

    def score_pose_pair(self, fp1, fp2, i=None):
        score = 0
        for r in fp1.feats:
            if r in fp2.feats:
                if i is not None:
                    score += fp1.feats[r][i]*fp2.feats[r][i]
                else:
                    score += np.dot(fp1.feats[r], fp2.feats[r])
        return score

    def score_all_pairs_of_poses(self):
        all_scores = {}
        for i1 in range(len(self.ligands)):
            for i2 in range(i1+1, len(self.ligands)):            
                poses1 = self.glides[self.ligands[i1]][self.struct].poses
                poses2 = self.glides[self.ligands[i2]][self.struct].poses

                # +-1 allows an extra (last) row/column to be used to score the crystallized poses
                edges = np.zeros((self.num_poses[self.ligands[i1]] + 1, self.num_poses[self.ligands[i2]] + 1)) 

                for p1 in range(-1,np.shape(edges)[0]-1):
                    for p2 in range(-1,np.shape(edges)[1]-1):
                        fp1 = poses1.get(p1, self.crystals[self.ligands[i1]]).fp
                        fp2 = poses2.get(p2, self.crystals[self.ligands[i2]]).fp
                        if fp1 == None or fp2 == None:
                            print self.ligands[i1], self.ligands[i2], self.struct, p1, p2
                        edges[p1][p2] = self.score_pose_pair(fp1, fp2)
            
                all_scores[(self.ligands[i1], self.ligands[i2])] = edges
    
        return all_scores
            
    def rank_poses_based_on_one_other_ligand(self):
        ranking = {}
        for l in self.ligands:
            ranking[l] = {}
            for (li, lj) in self.all_scores.keys():
                if l in (li, lj):
                    pair_scores = self.all_scores[(li, lj)] if l == li else self.all_scores[(li, lj)].T
                    # -1 excludes the crystal pose from comparison
                    scores = np.zeros( (3, np.shape(pair_scores)[0]) )
                    scores[1,:] = [np.max(pair_scores[p][:-1]) for p in range(np.shape(pair_scores)[0])]
                    scores[0,:] = [np.argmax(pair_scores[p][:-1]) for p in range(np.shape(pair_scores)[0])]
                    scores[2,:] = [self.get_rmsds(lj if l == li else li)[int(n)] for n in scores[0,:]]
                    ranking[l][lj if l == li else li] = scores
                
        return ranking

    def rank_poses_based_on_all_other_ligands(self):
        combined_ranking = {}
        for l in self.ligands:
            combined_ranking[l] = np.zeros(self.num_poses[l] + 1)
            for l2 in self.scores_from_one_ligand[l].keys():
                for p in range(len(combined_ranking[l])):
                    combined_ranking[l][p] += self.scores_from_one_ligand[l][l2][1][p]
    
        return combined_ranking

    def rmsd_of_neighbors(self, l, p):
        #final_pose = np.argmax(self.get_final_scores(l)[:-1])
        ave_rmsd = 0
        for l2 in self.scores_from_one_ligand[l].keys():
            ave_rmsd += self.scores_from_one_ligand[l][l2][2][p]/(len(self.ligands) - 1)
        return ave_rmsd

    def print_neighbors(self, l, p):
        for l2 in self.scores_from_one_ligand[l].keys():
            print l2, [self.scores_from_one_ligand[l][l2][i][p] for i in range(3)]

    def score_breakdown(self, l, p):
        #final_pose = np.argmax(self.get_final_scores(l)[:-1])
        fp_size = self.crystals[self.crystals.keys()[0]].fp.feats
        fp_size = len(fp_size[fp_size.keys()[0]])
        breakdown = [0 for i in range(fp_size)]
        for l2 in self.scores_from_one_ligand[l].keys():
            fp1 = self.glides[l][self.struct].poses[p].fp
            fp2 = self.glides[l2][self.struct].poses[self.scores_from_one_ligand[l][l2][0][p]].fp
            new = [self.score_pose_pair(fp1, fp2, i) for i in range(fp_size)]
            breakdown = [breakdown[i] + new[i] for i in range(fp_size)]
        return breakdown

    def get_rmsds(self,l):
        # returns all rmsds and then 0 for the crystal pose
        if l not in self.all_rmsds: 
            self.all_rmsds[l] = [self.glides[l][self.struct].poses[p].rmsd for p in range(self.num_poses[l])] + [0]
        return self.all_rmsds[l]

    def get_gscores(self,l):
        if l not in self.all_gscores:
            self.all_gscores[l] = [self.glides[l][self.struct].poses[p].gscore for p in range(self.num_poses[l])]
        return self.all_gscores[l]

    def get_all_scores(self,l1,l2):
        return self.all_scores[(l1,l2)] if (l1,l2) in self.all_scores else self.all_scores[(l2,l1)].T

    def get_final_scores(self, l):
        return self.final_scores_for_each_ligand[l]        

    def get_filtered_pose_pairs(self,l1,l2, rmsd_f, score_f, combo):
        all_scores = self.get_all_scores(l1,l2)
        rmsds1 = self.get_rmsds(l1)
        rmsds2 = self.get_rmsds(l2)
        out = []
        for p1 in range(len(rmsds1)-1):
            for p2 in range(len(rmsds2)-1):
                if combo(rmsd_f(rmsds1[p1],rmsds2[p2]), score_f(all_scores[p1][p2])):
                    out.append((p1, p2, rmsds1[p1], rmsds2[p2], all_scores[p1][p2]))
        return out

    def get_interactions(self, l1, l2, p1, p2):
        fp1 = self.glides[l1][self.struct].poses.get(p1, self.crystals[l1]).fp
        fp2 = self.glides[l2][self.struct].poses.get(p2, self.crystals[l2]).fp
        return self.score_pose_pair(fp1,fp2,get_details=True)

    def find_mismatched_interactions(self, int1, int2):
        '''
        returns list of mismatch, residue, interaction
        positive mismatches are stronger in int1
        negative mismatches are stronger in int2
        '''
        mismatch = []
        for (r, i) in int1:
            mismatch.append((int1[(r,i)] - int2.get((r,i), 0), r, i))
        for (r, i) in int2:
            if (r, i) not in int1:
                mismatch.append((-1*int2.get((r,i)), r, i))
        return mismatch


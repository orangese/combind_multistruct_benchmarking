import numpy as np

class Scores:
    def __init__(self, glides, crystals, ligands, structure, n, overlap=lambda x,y: x*y, weights={}):
        self.glides = glides
        self.crystals = crystals
        self.ligands = ligands
        self.struct = structure

        self.overlap = overlap
        self.w = weights
        
        self.num_poses = {l: min(n,len(glides[l][structure].poses)) for l in ligands}

        # the crystal poses are the last row/column of each matrix here
        self.all_scores = self.score_all_pairs_of_poses()
        self.scores_from_one_ligand = self.rank_poses_based_on_one_other_ligand()
        self.final_scores_for_each_ligand = self.rank_poses_based_on_all_other_ligands()

        self.all_rmsds = {}

    def score_pose_pair(self, fp1, fp2, get_details=False):
        score = 0
        interactions = {}
        for r in fp1.feats:
            if r in fp2.feats:
                for i in range(len(fp1.feats[r])):
                    i_score = self.overlap(fp1.feats[r][i]*self.w.get(i,1),fp2.feats[r][i]*self.w.get(i,1))
                    if get_details:
                        interactions[(r,i)] = i_score
                    score += i_score
        if get_details:
            return score, interactions
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
                    scores = [np.max(pair_scores[p][:-1]) for p in range(np.shape(pair_scores)[0])]
                    ranking[l][lj if l == li else li] = scores
                
        return ranking

    def rank_poses_based_on_all_other_ligands(self):
        combined_ranking = {}
        for l in self.ligands:
            combined_ranking[l] = np.zeros(self.num_poses[l] + 1)
            for l2 in self.scores_from_one_ligand[l].keys():
                for p in range(len(combined_ranking[l])):
                    combined_ranking[l][p] += self.scores_from_one_ligand[l][l2][p]
    
        return combined_ranking

    def get_rmsds(self,l):
        # returns all rmsds and then 0 for the crystal pose
        if l not in self.all_rmsds: 
            self.all_rmsds[l] = [self.glides[l][self.struct].poses[p].rmsd for p in range(self.num_poses[l])] + [0]
        return self.all_rmsds[l]

    def get_all_scores(self,l1,l2):
        return self.all_scores[(l1,l2)] if (l1,l2) in self.all_scores else self.all_scores[(l2,l1)].T

    def get_final_scores(self, l):
        return self.final_scores_for_each_ligand[l]        

    def get_filtered_pose_pairs(self,l1,l2, rmsd_f, score_f, combo):
        all_scores = self.get_all_scores(l1,l2)

        out = []
        for p1 in range(len(rmsds1)):
            for p2 in range(len(rmsds2)):
                if combo(rmsd_f(rmsds1[p1],rmsds2[p2]), score_f(all_scores[p1][p2])):
                    out.append(p1, p2, rmsds1[p1], rmsds2[p2], all_scores[p1][p2])
        return out

    def get_interactions(self, l1, l2, p1, p2):
        fp1 = self.glides[l1][self.struct].poses[p1]
        fp2 = self.glides[l2][self.struct].poses[p2]
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


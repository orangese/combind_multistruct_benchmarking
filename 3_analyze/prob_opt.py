import numpy as np

from pairs import LigPair

class PredictStructs:
    def __init__(self, docking_st, evidence, features, max_poses, T, reference):
        self.docking_st = docking_st
        self.ev = evidence
        self.features = features
        self.max_poses = max_poses
        self.T = float(T)
        self.reference = reference
        assert self.reference in ['DECOY', 'ALL', 'LTP']
        assert self.T >= 0

        self.ligand_partition_function_cache = {}
        self.log_likelihood_ratio_cache = {}
        self.lig_pairs = {}


    # Optimization algorithms
    def find_best_cluster(self, ligands, sampling=3, en_landscape=False):
        initial_cluster = {l:0 for l in ligands}
        return self._optimize_cluster(initial_cluster, sampling, en_landscape)

    def _optimize_cluster(self, pose_cluster, sampling, en_landscape):
        """
        Maximizes
            log P(L = l1 .. ln) = C + (\sum - glide / T) + (\sum log_odds)
        by the following method:

            1) Randomly select a ligand.
            2) Assume that the highest scoring pose for all other ligands is correct.
            3) Maximizing the score for the selected ligand.

        pose_cluster: dict of  ligand_name -> current pose number
        sampling: int sample sampling*len(init_cluster)**2 times
        """
        log_posteriors, rmsds = [], []

        time_since_update = 0
        while time_since_update < sampling * len(pose_cluster.keys())**2:
            time_since_update += 1

            query = np.random.choice(pose_cluster.keys())
            best_sc = self._partial_log_posterior(pose_cluster, query)
            best_p = pose_cluster[query]
            for p in range(self._num_poses(query)):
                pose_cluster[query] = p
                new_score = self._partial_log_posterior(pose_cluster, query)
                if new_score > best_sc:
                    best_sc, best_p = new_score, p
                    time_since_update = 0
            pose_cluster[query] = best_p

            if en_landscape:
                log_posteriors += [log_posterior(pose_cluster)]
                rmsds += [get_rmsd(pose_cluster)]

        return pose_cluster, (log_posteriors, rmsds)


    # Methods to score sets of ligands
    def log_posterior(self, pose_cluster):
        """
        Returns the log posterior for pose cluster.
        
        This should only be used for debugging purposes, as when optimizing we
        only need to consider terms that involve the ligand being considered.
        """
        log_prob = 0
        for ligname, pose in pose_cluster.items():
            log_prob += - self._get_gscore(ligname, pose) / self.T
        
        for i, (ligname1, pose1) in enumerate(pose_cluster.items()):
            for ligname2, pose2 in pose_cluster.items()[i+1:]:
                log_prob += self._log_likelihood_ratio_pair(pose_cluster, ligname1, ligname2)
        return energy

    def _partial_log_posterior(self, pose_cluster, query):
        """
        Computes the partial contribution of ligand 'query' in 'pose' to the total log prob.
        This is simply - GlideScore(l_q)  / T + \sum log_odds.
        """
        log_odds = sum(self._log_likelihood_ratio_pair(pose_cluster, query, ligname)
                       for ligname in pose_cluster if ligname != query)
        log_prior = - self._get_gscore(query, pose_cluster[query]) / self.T
        return log_odds + log_prior

    def _log_likelihood_ratio_pair(self, pose_cluster, ligname1, ligname2):
        """
        Computes the pairwise log likelihood ratio for the poses

                poses_cluster[ligname1] and pose_clusters[ligname2]
        
        as the the sum of the log likelihoods of each of the features.

        This function supports 3 options for the "denominator" or "reference"
        distribution that is specified by self.reference.
        """
        pair_key = ((ligname1, ligname2, pose_cluster[ligname1], pose_cluster[ligname2])
                    if ligname1 < ligname2 else
                    (ligname2, ligname1, pose_cluster[ligname2], pose_cluster[ligname1]))
        
        if pair_key not in self.log_likelihood_ratio_cache:
            log_likelihood = 0
            for k, k_def in self.features.items():

                k_val = self._get_feature(k, ligname1, ligname2,
                                          pose_cluster[ligname1], pose_cluster[ligname2])
                
                if k_val is None: continue
                assert k_val <= 1 and k_val >= 0, "{} {}".format(k, k_val)

                p_x_native  = self.ev.evaluate(k, k_val, 1)
                
                # Choose between 3 potential options for reference distribution
                if self.reference == 'DECOY':
                    p_x= self.ev.evaluate(k, k_val, 0)
                elif self.reference == 'ALL':
                    p_x = self.ev.evaluate(k, k_val, -1)   
                elif self.reference == 'LTP':
                    prior = (  self._get_prior(ligname1, pose_cluster[ligname1])
                             * self._get_prior(ligname2, pose_cluster[ligname2]))
                    p_x = (  self.ev.evaluate(k, k_val, 1) * prior
                           + self.ev.evaluate(k, k_val, 0) * (1 - prior))            

                log_likelihood += np.log(p_x_native) - np.log(p_x)
            self.log_likelihood_ratio_cache[pair_key] = log_likelihood
        return self.log_likelihood_ratio_cache[pair_key]

    def _get_feature(self, fname, l1, l2, p1, p2):
        if l1 > l2:
            l1, l2 = l2, l1
            p1, p2 = p2, p1

        if (l1, l2) not in self.lig_pairs:

            self.lig_pairs[(l1, l2)] = LigPair(self.docking_st.ligands[l1], self.docking_st.ligands[l2],
                                               self.features, self.docking_st.mcss)
        return self.lig_pairs[(l1, l2)].get_feature(fname, p1, p2)


    # Methods related to physics scores
    def _get_prior(self, ligname, pose):
        return ( np.exp(- self._get_gscore(ligname, pose) / self.T)
                / self._get_ligand_partition_function(ligname))

    def _get_ligand_partition_function(self, ligname):
        if ligname not in self.ligand_partition_function_cache:
            self.ligand_partition_function_cache[ligname] = sum(
                np.exp(- self._get_gscore(ligname, pose) / self.T)
                for pose in range(self._num_poses(ligname)))
        return self.ligand_partition_function_cache[ligname]

    def _get_gscore(self, ligand, pose):
        return self.docking_st.ligands[ligand].poses[pose].gscore


    # Methods to manage poses
    def get_poses(self, cluster):
        return {l:self.docking_st.ligands[l].poses[p] for l,p in cluster.items()}

    def get_rmsd(self, cluster):
        return np.mean([self.docking_st.ligands[l].poses[p].rmsd for l,p in cluster.items()])

    def _num_poses(self, ligname):
        return min(self.max_poses, len(self.docking_st.ligands[ligname].poses))


    ############################################################################

    # def x(self, pose_cluster, k, k_def, sample_l=None, lig_order=None):
    #     if lig_order is None:
    #         lig_order = pose_cluster.keys()

    #     log_p_l_given_x = np.zeros( (len(lig_order), len(lig_order)) )
    #     for i, l1 in enumerate(lig_order):
    #         for j, l2 in list(enumerate(lig_order))[i+1:]:
    #             if sample_l is not None and sample_l != l1 and sample_l != l2:
    #                 continue
    #             #pp = self.ligset.get_pose_pair(l1, pose_cluster[l1], l2, pose_cluster[l2])
    #             lp = self.ligset.get_lig_pair(l1, l2)
    #             k_val = lp.get_feature(k, pose_cluster[l1],pose_cluster[l2])
    #             if k_val is None: continue
    #             assert k_val <= 1 and k_val >= 0, '{} {}'.format(fname, kval)
    #             x[i,j] = k_val
         
    #             prior = np.exp(self.ligset.log_prior(l1,pose_cluster[l1])+self.ligset.log_prior(l2,pose_cluster[l2]))
       
    #             p_x_native = self.ligset.ev.evaluate(k, k_val, 1)
    #             p_x_nnative = self.ligset.ev.evaluate(k, k_val, 0)

    #             # option A:
    #             p_x = p_x_native*prior + p_x_nnative*(1-prior)

    #             # option B:
    #             #p_x = self.ligset.ev.evaluate(k, k_val, -1)               

    #             log_p_l_given_x[i,j] = np.log(p_x_native/p_x)
    #             if p_x == 0:
    #                 print l1,l2,k, k_val
    #                 print p_x_native, p_x

    #     return log_p_l_given_x




    # def joint_posterior(self, pose_cluster, sample_l=None, verbose=False, en_landscape=False):
    #     # pose cluster = ligand_name: integer_pose_index
        
    #     log_prob = sum([self.ligset.log_prior(l,p) for l,p in pose_cluster.items()])
    #     for k,k_def in self.ligset.features.items():
    #         x_k, log_p_l_given_x_k = self.x(pose_cluster, k, k_def, sample_l=sample_l)
    #         log_prob += np.sum(log_p_l_given_x_k)

    #     if en_landscape:
    #         return log_prob, self.ligset.get_rmsd(pose_cluster)
    #     return log_prob, None

    # def max_posterior(self, ligands, verbose=False, sampling=3, en_landscape=False):
    #     initial_cluster = {l:0 for l in ligands}
    #     max_sc, best_cluster, all_clusters = self.optimize(initial_cluster,sampling=sampling, en_landscape=en_landscape)

    #     if verbose:
    #         print 'cluster -1, score {}'.format(max_sc)#, max_sc#, rmsd:', self.joint_posterior(best_cluster, en_landscape=True)

    #     # run 10 more times starting randomly
    #     for i in range(15):
    #         rand_cluster = {}
    #         for l in ligands:
    #             rand_cluster[l] = np.random.randint(self.ligset.num_poses[l])

    #         new_sc, new_cluster, clusters = self.optimize(rand_cluster,sampling=sampling, en_landscape=en_landscape)
    #         all_clusters.extend(clusters)

    #         if new_sc > max_sc:
    #             max_sc, best_cluster = new_sc, new_cluster
    #         if verbose:
    #             print 'cluster {}, score {}'.format(i,new_sc)# rmsd:'.format(i), self.joint_posterior(pose_cluster, en_landscape=True)

    #     return best_cluster, all_clusters



    # def optimize(self, init_cluster, sampling=3, en_landscape=False):
    #     # maximize objective by randomly sampling ligands
    #     time_since_update = 0
    #     pose_cluster = init_cluster

    #     all_clusters = []

    #     while time_since_update < sampling*(len(init_cluster.keys()))**2:
    #         time_since_update += 1

    #         rand_lig = np.random.choice(init_cluster.keys())
    #         sample_l = rand_lig
    #         if en_landscape: sample_l = None

    #         max_sc, rmsd = self.joint_posterior(pose_cluster, sample_l=sample_l, en_landscape=en_landscape)
    #         best_p = pose_cluster[rand_lig]
    #         if en_landscape: all_clusters.append((max_sc, rmsd))
            
    #         for p in range(self.ligset.num_poses[rand_lig]):
    #             pose_cluster[rand_lig] = p
    #             new_score, new_rmsd = self.joint_posterior(pose_cluster, sample_l=sample_l, en_landscape=en_landscape)
    #             if new_score > max_sc:
    #                 max_sc, best_p = new_score, p
    #                 time_since_update = 0
    #             if en_landscape: all_clusters.append((new_score, new_rmsd))

    #         pose_cluster[rand_lig] = best_p

    #     return self.joint_posterior(pose_cluster)[0], pose_cluster, all_clusters
            
    # def score_query(self, query, optimal_cluster):
    #     eval_cluster = {l:p for l,p in optimal_cluster.items()}
    #     scores = []
    #     for p in range(self.ligset.n(query)):
    #         eval_cluster[query] = p
    #         scores.append(self.joint_posterior(eval_cluster, sample_l=query)[0])
    #     return scores




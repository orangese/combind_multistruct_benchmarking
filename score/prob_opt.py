"""
Core optimization code.
"""

import numpy as np
from score.pairs import LigPair

class PredictStructs:
    """
    ligands ({name: containers.Ligand, })
    mcss (mcss.MCSSController): All applicable RMSDs should be loaded.
    stats ({feature: {'native': score.DensityEstimate, 'reference': score.DensityEstimate}})
    features ([str, ]): Features to use when computing similarity scores.
    max_poses (int): Maximum number of poses to consider.
    alpha (float): Factor to which to weight the glide scores.
    """
    def __init__(self, ligands, mcss, stats, features, max_poses, alpha):
        self.ligands = ligands
        self.mcss = mcss
        self.stats = stats
        self.features = features
        self.max_poses = max_poses
        self.alpha = float(alpha)

        self.log_likelihood_ratio_cache = {}
        self.lig_pairs = {}

        self._validate()

    def _validate(self):
        if 'mcss' in self.features:
            assert self.mcss is not None
            for lig1 in self.ligands:
                for lig2 in self.ligands:
                    if lig1 == lig2: continue
                    try:
                        self.mcss.get_rmsd(lig1, lig2, 0, 0)
                    except KeyError:
                        assert False
        
        for feature in self.features:
            assert feature in self.stats['native']
            assert feature in self.stats['reference']
        
        assert type(self.ligands) == dict
        assert self.alpha >= 0
        assert self.max_poses > 0

    def max_posterior(self, max_iterations = 1e6, restart = 50):
        """
        Computes the pose cluster maximizing the posterior likelihood.

        Note that this function is not guarenteed to find the global maximum
        since the search space is not convex. This is handled by restarting from
        multiple initial configurations and returning the best scoring configuration
        found in any optimization.

        max_iterations (int): Maximum number of iterations to attempt before exiting.
        restart (int): Number of times to run the optimization
        """
        self._validate()
        best_score = -float('inf')
        for i in range(restart):
            cluster = ({lig: 0 for lig in self.ligands}
                       if i == 0 else
                       {lig: np.random.randint(self._num_poses(lig)) for lig in self.ligands})

            score, cluster = self._optimize_cluster(cluster, max_iterations)
            if score > best_score:
                best_score = score
                best_cluster = cluster
           
            print(cluster)
            print('cluster {}, score {}'.format(i, score))

        return best_cluster

    def _optimize_cluster(self, pose_cluster, max_iterations):
        """
        Maximizes
            log P(L = l1 .. ln) = C + (sum - glide / T) + (sum log_odds)
        by the following method:

            1) Randomly select a ligand.
            2) Assume that the highest scoring pose for all other ligands is correct.
            3) Maximizing the score for the selected ligand.

        pose_cluster ({ligand_name: current pose number, })
        max_iterations (int)
        """
        iterations = 0
        update = True
        while update:
            if iterations > max_iterations:
                assert False, ('Warning: Reached maximum iterations'
                               '({}) without convergence.'.format(max_iterations))
            iterations += 1
            
            update = False
            for query in np.random.permutation(list(pose_cluster.keys())):
                best_score = self._partial_log_posterior(pose_cluster, query)
                best_pose  = pose_cluster[query]
                for pose in range(self._num_poses(query)):
                    pose_cluster[query] = pose
                    score = self._partial_log_posterior(pose_cluster, query)
                    if score > best_score:
                        best_score, best_pose = score, pose
                        update = True
                pose_cluster[query] = best_pose

        return self.log_posterior(pose_cluster), pose_cluster

    # Methods to score sets of ligands
    def log_posterior(self, pose_cluster):
        """
        Returns the log posterior for pose cluster.
        
        This should only be used for debugging purposes, as when optimizing we
        only need to consider terms that involve the ligand being considered.
        """
        log_prob = 0
        for ligname, pose in pose_cluster.items():
            log_prob -= self._get_gscore(ligname, pose) * self.alpha
        
        ligands = list(pose_cluster.keys())
        for i, ligname1 in enumerate(ligands):
            for ligname2 in ligands[i+1:]:
                log_prob += self._log_likelihood_ratio_pair(pose_cluster, ligname1, ligname2)
        return log_prob

    def _partial_log_posterior(self, pose_cluster, query):
        """
        Computes the partial contribution of ligand 'query' in 'pose' to the total log prob.
        This is simply - GlideScore(l_q)  / T + sum log_odds.
        """
        log_odds = sum(self._log_likelihood_ratio_pair(pose_cluster, query, ligname)
                       for ligname in pose_cluster if ligname != query)
        log_prior = - self._get_gscore(query, pose_cluster[query]) * self.alpha
        return log_odds + log_prior

    def _log_likelihood_ratio_pair(self, pose_cluster, ligname1, ligname2):
        """
        Computes the pairwise log likelihood ratio for the poses

                poses_cluster[ligname1] and pose_clusters[ligname2]
        
        as the the sum of the log likelihoods of each of the k_list.
        """
        pair_key = ((ligname1, ligname2, pose_cluster[ligname1], pose_cluster[ligname2])
                    if ligname1 < ligname2 else
                    (ligname2, ligname1, pose_cluster[ligname2], pose_cluster[ligname1]))
        
        if pair_key not in self.log_likelihood_ratio_cache:
            log_likelihood = 0
            for feature in self.features:
                _, p_x_native, p_x = self._likelihoods_for_pair_and_single_feature(feature,
                                                                                   pose_cluster,
                                                                                   ligname1,
                                                                                   ligname2)
                # If either of these are 0, we will get infs or nans
                # using a gKDE this should not happen, but if we ever switch
                # to a compact kernel, this could be hard to debug
                assert p_x_native != 0.0
                assert p_x != 0.0

                log_likelihood += np.log(p_x_native) - np.log(p_x)
            self.log_likelihood_ratio_cache[pair_key] = log_likelihood
        return self.log_likelihood_ratio_cache[pair_key]

    def _likelihoods_for_pair_and_single_feature(self, feature, pose_cluster, ligname1, ligname2):
        """
        Returns the feature value 'x' and its likelihoods P(x|l) and P(x)
        for feature 'fname' and poses pose_cluster[ligname1] and pose_cluster[ligname2].
        """
        x = self._get_feature(feature, ligname1, ligname2,
                              pose_cluster[ligname1], pose_cluster[ligname2])
        
        if x is None: return 0.0, 1.0, 1.0

        p_x_native  = self.stats['native'][feature](x)
        p_x = self.stats['reference'][feature](x)

        return x, p_x_native, p_x

    def _get_feature(self, feature, ligname1, ligname2, pose1, pose2):
        # Maintain the convention that ligname1 < ligname2, so that
        # we don't put duplicate entries in self.lig_pairs.
        if ligname1 > ligname2:
            ligname1, ligname2 = ligname2, ligname1
            pose1, pose2 = pose2, pose1

        # Constructing LigPairs is expensive, so we cache them.
        if (ligname1, ligname2) not in self.lig_pairs:
            self.lig_pairs[(ligname1, ligname2)] = LigPair(self.ligands[ligname1],
                                                           self.ligands[ligname2],
                                                           self.features, self.mcss,
                                                           self.max_poses)

        return self.lig_pairs[(ligname1, ligname2)].get_feature(feature, pose1, pose2)

    def _get_gscore(self, ligname, pose):
        return self.ligands[ligname].poses[pose].gscore

    def _num_poses(self, ligname):
        return min(self.max_poses, len(self.ligands[ligname].poses))

    ####################################################################################
    # Methods important for figure making, but not execution.
    def likelihood_and_feature_matrix(self, pose_cluster, k, lig_order):
        """
        Returns the feature values and likelihood ratios, P(X|l)/P(X)
        for feature 'fname' for poses in 'pose_cluster'
        as len(pose_cluster) x len(pose_cluster) numpy arrays.
        """
        x                    = np.zeros((len(pose_cluster), len(pose_cluster)))
        log_likelihood_ratio = np.zeros((len(pose_cluster), len(pose_cluster)))

        for i, ligname1 in enumerate(lig_order):
            for j, ligname2 in enumerate(lig_order):
                if j <= i: continue
                x_k, p_x_native, p_x = self._likelihoods_for_pair_and_single_feature(k, pose_cluster,
                                                                                   ligname1, ligname2)
                x[i, j] = x_k
                log_likelihood_ratio[i, j] = np.log(p_x_native) - np.log(p_x)
        return x, log_likelihood_ratio
    
    def get_poses(self, cluster):
        return {l:self.docking[l].poses[p] for l,p in cluster.items()}

    def get_rmsd(self, cluster):
        tmp = [self.docking[l].poses[p].rmsd for l,p in cluster.items()]
        return np.mean([r for r in tmp if r is not None])

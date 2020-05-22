"""
Core optimization code.
"""

import numpy as np
from score.pairs import LigPair

class PredictStructs:
    """
    stats ({feature: {'native': score.DensityEstimate, 'reference': score.DensityEstimate}})
    features ([str, ]): Features to use when computing similarity scores.
    max_poses (int): Maximum number of poses to consider.
    alpha (float): Factor to which to weight the glide scores.
    """
    def __init__(self, stats, features, max_poses, alpha):
        self.stats = stats
        self.features = features
        self.max_poses = max_poses
        self.alpha = float(alpha)

        self.ligands = {}
        self.mcss = None
        self.log_likelihood_ratio_cache = {}
        self.lig_pairs = {}

    def set_ligands(self, ligands, mcss):
        self.ligands = ligands
        self.mcss = mcss
        self.log_likelihood_ratio_cache = {}
        self.lig_pairs = {}

    def max_posterior(self, max_iterations, restart):
        """
        Computes the pose cluster maximizing the posterior likelihood.

        Note that this function is not guarenteed to find the global maximum
        since the search space is not convex. This is handled by restarting from
        multiple initial configurations and returning the best scoring configuration
        found in any optimization.

        max_iterations (int): Maximum number of iterations to attempt before exiting.
        restart (int): Number of times to run the optimization
        """
        best_score = -float('inf')
        for i in range(restart):
            if i == 0:
                cluster = {lig: 0 for lig in self.ligands}
            else:
                cluster = {lig: np.random.randint(self._num_poses(lig))
                           for lig in self.ligands}

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

        # Core coordinate ascent logic.
        for _ in range(max_iterations):
            update = False
            for query in np.random.permutation(list(pose_cluster.keys())):
                best_pose = self._best_pose(pose_cluster, query)
                if best_pose != pose_cluster[query]:
                    update = True
                    pose_cluster[query] = best_pose
            if not update:
                break
        return self.log_posterior(pose_cluster), pose_cluster


    def _best_pose(self, pose_cluster, query):
        # Copy the pose cluster so that we don't have to worry
        # about mutating it.
        pose_cluster = {k:v for k, v in pose_cluster.items()}
        best_score = -float('inf')
        best_pose  = -1
        for pose in range(self._num_poses(query)):
            pose_cluster[query] = pose
            score = self._partial_log_posterior(pose_cluster, query)
            if score > best_score:
                best_score, best_pose = score, pose
        return best_pose

    # Methods to score sets of ligands
    def log_posterior(self, pose_cluster):
        """
        Returns the log posterior for pose cluster.
        
        This should only be used for debugging purposes, as when optimizing we
        only need to consider terms that involve the ligand being considered.
        """
        log_prob = 0
        for ligname, pose in pose_cluster.items():
            log_prob -= (self._get_physics_score(ligname, pose) * self.alpha
                         * self._effective_number(pose_cluster, ligname))
        
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
        log_odds = 0
        for ligname in pose_cluster:
            if ligname == query: continue
            log_odds += self._log_likelihood_ratio_pair(pose_cluster, query, ligname)
        
        log_prior = -self._get_physics_score(query, pose_cluster[query])
        log_prior *= self.alpha * self._effective_number(pose_cluster, query)
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
                # to a compact kernel, this could be hard to debug.
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

    # Helpers
    def _effective_number(self, pose_cluster, query):
        return sum(1 for ligand, pose in pose_cluster.items()
                   if ligand != query)

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

    def _get_physics_score(self, ligname, pose):
        return self.ligands[ligname].poses[pose].gscore

    def _num_poses(self, ligname):
        return min(self.max_poses, len(self.ligands[ligname].poses))

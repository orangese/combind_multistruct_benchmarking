import numpy as np
import random
from pairs import LigPair

class PredictStructs:
    def __init__(self, docking_st, evidence, features, max_poses, T, reference='ALL', normalize_fp=True):
        self.docking_st = docking_st
        self.ev = evidence
        self.features = features
        self.max_poses = max_poses
        self.T = float(T)
        self.reference = reference
        self.normalize_fp = normalize_fp
        assert self.reference in ['DECOY', 'ALL', 'LTP', 'OLD']
        assert self.T >= 0
        assert self.normalize_fp == evidence.normalize_fp

        self.ligand_partition_function_cache = {}
        self.log_likelihood_ratio_cache = {}
        self.lig_pairs = {}

    def max_posterior(self, ligands, verbose=False, sampling=3, mode='MAX', restart=15, en_landscape=False):

        opt = self._optimize_cluster
        if mode == 'ANNEAL':
            opt = self._anneal_cluster

        initial_cluster = {l:0 for l in ligands}
        max_sc, best_cluster, all_scores, all_rmsds = opt(initial_cluster,sampling=sampling, en_landscape=en_landscape)

        if verbose:
            print 'cluster -1, score {}'.format(max_sc)

        for i in range(restart):
            rand_cluster = {}
            for l in ligands:
                rand_cluster[l] = np.random.randint(self._num_poses(l))

            new_sc, new_cluster, scores, rmsds = opt(rand_cluster,sampling=sampling, en_landscape=en_landscape)
            all_scores.extend(scores)
            all_rmsds.extend(rmsds)

            if new_sc > max_sc:
                max_sc, best_cluster = new_sc, new_cluster
            if verbose:
                print 'cluster {}, score {}'.format(i,new_sc)

        return best_cluster, all_scores, all_rmsds

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
                log_posteriors += [self.log_posterior(pose_cluster)]
                rmsds += [self.get_rmsd(pose_cluster)]

        return self.log_posterior(pose_cluster), pose_cluster, log_posteriors, rmsds

    def _anneal_cluster(self, pose_cluster, sampling, en_landscape):
        """
        Optimize pose cluster using simulated annealing.

        This method has been tested, but I haven't put any effort
        into tuning the cooling schedule. It isn't immediately way better.

        Keeps track of the score of the current pose cluster by
        summing the difference in partial score at each switch.
        Should verify that this doesn't lead to numerical issues,
        but I really doubt it.
        """
        
        log_likelihood = 0 # confusing to call this energy because we are maximizing it
        best_log_likelihood = 0
        best_cluster = {k:v for k, v in pose_cluster.items()}
        
        # Linear cooling schedule
        T_START = 2
        COOLING = 0.5
        CHAIN_LENGTH = sampling*self.max_poses*len(pose_cluster)
        T = lambda i: float(T_START - COOLING*int(i/CHAIN_LENGTH))

        for i in range(CHAIN_LENGTH*int(T_START/COOLING) + 1):

            query = np.random.choice(pose_cluster.keys())
            current_pose  = pose_cluster[query]
            proposed_pose = random.randint(0, self._num_poses(query)-1)

            current_partial  = self._partial_log_posterior(pose_cluster, query)
            pose_cluster[query] = proposed_pose
            proposed_partial = self._partial_log_posterior(pose_cluster, query)
        
            if (proposed_partial > current_partial
                or T(i) > 0 and random.random() < np.exp((proposed_partial - current_partial) / T(i))):
                pose_cluster[query] = proposed_pose
                log_likelihood += proposed_partial - current_partial
            else:
                pose_cluster[query] = current_pose
                
            if log_likelihood > best_log_likelihood:
                best_log_likelihood = log_likelihood
                best_cluster = {k:v for k, v in pose_cluster.items()}

        return self.log_posterior(pose_cluster), best_cluster, [], []

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
        return log_prob

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
                k_val, p_x_native, p_x = self._likelihoods_for_pair_and_single_feature(k, pose_cluster,
                                                                                   ligname1, ligname2)
                x[i, j] = k_val
                log_likelihood_ratio[i, j] = np.log(p_x_native) - np.log(p_x)
        return x, log_likelihood_ratio

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
            for k in self.features:
                _, p_x_native, p_x = self._likelihoods_for_pair_and_single_feature(k, pose_cluster,
                                                                                   ligname1, ligname2)
                log_likelihood += np.log(p_x_native) - np.log(p_x)
            self.log_likelihood_ratio_cache[pair_key] = log_likelihood
        return self.log_likelihood_ratio_cache[pair_key]


    def _likelihoods_for_pair_and_single_feature(self, k, pose_cluster, ligname1, ligname2):
        """
        Returns the feature value 'x' and its likelihoods P(x|l) and P(x)
        for feature 'fname' and poses pose_cluster[ligname1] and pose_cluster[ligname2].
        """
        k_val = self._get_feature(k, ligname1, ligname2,
                                  pose_cluster[ligname1], pose_cluster[ligname2])
        
        if k_val is None: return 0, 1, 1
        if self.normalize_fp:
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
        elif self.reference == 'OLD':
            p_x_native = np.exp(k_val) if k != 'mcss' else 1.0
            p_x = 1.0

        return k_val, p_x_native, p_x

    def score_query(self, query, optimal_cluster):
         print 'fixme'
#        eval_cluster = {l:p for l,p in optimal_cluster.items()}
#        scores = []
#        for p in range(self.ligset.n(query)):
#            eval_cluster[query] = p
#            scores.append(self.joint_posterior(eval_cluster, sample_l=query)[0])
#        return scores

    def _get_feature(self, fname, l1, l2, p1, p2):
        if l1 > l2:
            l1, l2 = l2, l1
            p1, p2 = p2, p1

        if (l1, l2) not in self.lig_pairs:

            self.lig_pairs[(l1, l2)] = LigPair(self.docking_st.ligands[l1], self.docking_st.ligands[l2],
                                               self.features, self.docking_st.mcss, self.max_poses, self.normalize_fp)
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
        tmp = [self.docking_st.ligands[l].poses[p].rmsd for l,p in cluster.items()]
        return np.mean([r for r in tmp if r is not None])

    def _num_poses(self, ligname):
        return min(self.max_poses, len(self.docking_st.ligands[ligname].poses))

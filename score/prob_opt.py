"""
Core optimization code.
"""

import numpy as np
from score.pairs import LigPair
from matplotlib import colors
from matplotlib.gridspec import GridSpec
from matplotlib.colors import LinearSegmentedColormap
import matplotlib.pyplot as plt
from shared_paths import feature_defs

class PredictStructs:
    """
    ligands ({name: containers.Ligand, })
    mcss (mcss.MCSSController): All applicable RMSDs should be loaded.
    stats ({feature: {'native': score.DensityEstimate, 'reference': score.DensityEstimate}})
    features ([str, ]): Features to use when computing similarity scores.
    max_poses (int): Maximum number of poses to consider.
    alpha (float): Factor to which to weight the glide scores.
    overlap (float): Overlap metric. Must match stats.
    xtal_boost (float): Multiplier for pair scores involving xtal poses.
    physics_score ('gscore' or 'emodel'): which physics score to use.
    all_wrong (bool): Whether to allow possibility that all poses are wrong
        and exclude the ligand from the optimization.
    """
    def __init__(self, ligands, mcss, stats,
                 features, max_poses, alpha, overlap='maxoverlap',
                 xtal_boost=1.0, physics_score='gscore', all_wrong=False):
        self.ligands = ligands
        self.mcss = mcss
        self.stats = stats
        self.features = features
        self.overlap = overlap
        self.max_poses = max_poses
        self.alpha = float(alpha)
        self.xtal_boost = xtal_boost
        self.physics_score = physics_score
        self.all_wrong = -1 if all_wrong else 0

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

    def max_posterior(self, max_iterations=1000000, restart=500):
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
                best_pose = self._best_pose(pose_cluster, query, self.all_wrong)
                if best_pose != pose_cluster[query]:
                    update = True
                    pose_cluster[query] = best_pose
            if not update:
                break

        # Assign all poses in optimal cluster to best real pose.
        output = {}
        for query, pose in pose_cluster.items():
            if pose != -1:
                output[query] = pose
            else:
                output[query] =self._best_pose(pose_cluster, query, 0)
        return self.log_posterior(output), output


    def _best_pose(self, pose_cluster, query, all_wrong):
        # Copy the pose cluster so that we don't have to worry
        # about mutating it.
        pose_cluster = {k:v for k, v in pose_cluster.items()}
        best_score = -float('inf')
        best_pose  = -1
        for pose in range(all_wrong, self._num_poses(query)):
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
            log_prob -= self._get_physics_score(ligname, pose) * self.alpha
        
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
            lr = self._log_likelihood_ratio_pair(pose_cluster, query, ligname)
            
            lr *= self.xtal_boost if 'crystal' in ligname else 1.0
            log_odds += lr
        
        log_prior = -self._get_physics_score(query, pose_cluster[query])
        log_prior *= self.alpha * self._effective_number(pose_cluster, query)
        return log_odds + log_prior

    def _log_likelihood_ratio_pair(self, pose_cluster, ligname1, ligname2):
        """
        Computes the pairwise log likelihood ratio for the poses

                poses_cluster[ligname1] and pose_clusters[ligname2]
        
        as the the sum of the log likelihoods of each of the k_list.
        """
        if pose_cluster[ligname1] == -1 or pose_cluster[ligname2] == -1:
            return 0

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
                   if pose != -1 and ligand != query)

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
                                                           self.max_poses,
                                                           self.overlap)

        return self.lig_pairs[(ligname1, ligname2)].get_feature(feature, pose1, pose2)

    def _get_physics_score(self, ligname, pose):
        if self.physics_score == 'gscore':
            if pose == -1:
                return -8.0
            return self.ligands[ligname].poses[pose].gscore
        else:
            if pose == -1:
                assert False
                return -80.0 # Not sure if this is right?
            return self.ligands[ligname].poses[pose].emodel

    def _num_poses(self, ligname):
        return min(self.max_poses, len(self.ligands[ligname].poses))


class PredictStructsFigures(PredictStructs):
    """
    Methods of this class are used for figure making. They are seperated
    from the PredictStructs to prevent obscuring the core optimization logic.
    """
    def likelihood_and_feature_matrix(self, pose_cluster, k, lig_order):
        """
        Returns the feature values and likelihood ratios, P(X|l)/P(X)
        for feature 'fname' for poses in 'pose_cluster'
        as len(pose_cluster) x len(pose_cluster) numpy arrays.
        """
        x                    = np.zeros((len(lig_order), len(lig_order)))
        log_likelihood_ratio = np.zeros((len(lig_order), len(lig_order)))

        for i, ligname1 in enumerate(lig_order):
            for j, ligname2 in enumerate(lig_order):
                if j <= i: continue
                x_k, p_x_native, p_x = self._likelihoods_for_pair_and_single_feature(k, pose_cluster,
                                                                                   ligname1, ligname2)
                x[i, j] = x_k
                log_likelihood_ratio[i, j] = np.log(p_x_native) - np.log(p_x)
        return x, log_likelihood_ratio

    def interactions(self, pose_clusters, k):
        interactions = set()
        for pose_cluster in pose_clusters:
            _interactions = set([interaction
                                 for ligand, pose in pose_cluster.items()
                                 for interaction in self.ligands[ligand].poses[pose].fp
                                 if (interaction[0] in feature_defs[k])])
            interactions = interactions.union(_interactions)
        interactions = sorted(interactions)
        X = 0
        for pose_cluster in pose_clusters:
            X += self.interaction_matrix(pose_cluster, interactions, list(pose_cluster.keys()))

        # Remove interactions that are never seen
        interactions = [interaction
                        for interaction, present in zip(interactions, X.sum(axis = 1) > 0)
                        if present]
        X = X[X.sum(axis = 1) > 0]

        # Sort by interacton frequency
        idx = np.argsort(X.sum(axis = 1))[::-1]
        return [interactions[i] for i in idx]

    def interaction_matrix(self, pose_cluster, interactions, lig_order):
        X = []
        for ligand in lig_order:
            X += [[]]
            pose = self.ligands[ligand].poses[pose_cluster[ligand]]
            for interaction in interactions:
                X[-1] += [pose.fp[interaction] if interaction in pose.fp else 0]
        return np.array(X).T

    def gel_plot_individual(self, cluster, k, lig_order, interactions=None,
                            divide=2, resname_size=14, ax=None, pretty=True):
        if interactions is None:
            interactions = self.interactions([cluster], k)
        X = self.interaction_matrix(cluster, interactions, lig_order)
        labels = [self._format_int(interaction) for interaction in interactions]

        # Y is X + empty cells seperating boxes.
        Y = np.zeros((X.shape[0]*divide, X.shape[1]*divide))
        Y[::divide, ::divide] = X
        if ax is None:
            figure = plt.figure(figsize =  (9, 5))
            gs = GridSpec(1, 2, figure = figure, width_ratios = [4, 5])
            ax = plt.subplot(gs[1])
        ax.imshow(Y, cmap='binary', aspect = 'auto', vmin=0, vmax=1.0)
        #ax.set(adjustable='box', aspect='equal')

        if pretty:
            self._pretty(lig_order, labels, divide, resname_size)

    
    def gel_plot(self, combind_cluster, glide_cluster, k, lig_order,
                 divide=2, resname_size=14):
        interactions = self.interactions([combind_cluster, glide_cluster], k)
        combind = self.interaction_matrix(combind_cluster,interactions, lig_order)
        glide = self.interaction_matrix(glide_cluster, interactions, lig_order)

        seen = combind + glide
        X = combind - glide

        labels = [self._format_int(interaction) for interaction in interactions]

        # Differentiate values that are the same and
        # positive from those that are zero.
        X[seen == 0] = float('inf')
    
        # Y is X + empty cells seperating boxes.
        Y = np.zeros((X.shape[0]*divide, X.shape[1]*divide))
        Y.fill(float('inf')) # dividers should be white.
        Y[::divide, ::divide] = X

    
        rkb = [(1, 0, 0.5), (0, 0, 0), (0.5, 1, 0)]
        cm = LinearSegmentedColormap.from_list('rkb', rkb, N=256, gamma=1.0)
        figure = plt.figure(figsize =  (9, 5))
        gs = GridSpec(1, 2, figure = figure, width_ratios = [4, 5])
        ax = plt.subplot(gs[1])
        ax.imshow(Y, cm, vmin=-0.5, vmax=0.5)
        ax.set(adjustable='box', aspect='equal')
        self._pretty(lig_order, labels, divide, resname_size)
        print('ComBind - Glide: {}'.format(X.sum()))

    def _pretty(self, lig_order, resnames, divide, resname_size):
        for spine in plt.gca().spines.values():
            spine.set_visible(False)
        for axis in ['x', 'y']:
            plt.tick_params(axis=axis, which='both',bottom=False,
                            top=False, left=False, right=False)
        plt.yticks(range(0, divide*len(resnames), divide), resnames,
                   size = resname_size, fontname = 'monospace')
        plt.xticks(range(0, divide*len(lig_order), divide), lig_order,
                   rotation = 'vertical', size = resname_size, fontname = 'monospace')

    def _format_int(self, interaction):
        three_to_one = {
            'TYR': 'Y',
            'VAL': 'V',
            'ARG': 'R',
            'THR': 'T',
            'GLU': 'E',
            'SER': 'S',
            'PHE': 'F',
            'ALA': 'A',
            'MET': 'M',
            'ILE': 'I',
            'ASP': 'D',
            'GLN': 'Q',
            'ASN': 'N',
            'GLY': 'G',
            'PRO': 'P',
            'TRP': 'W',
            'LEU': 'L',
            'LYS': 'K',
            'CYS': 'C',
            'HIS': 'H'
        }
        name = interaction[1].split('(')[1][:-1]
        num = interaction[1].split('(')[0].split(':')[1]
        if name in three_to_one:
            name = three_to_one[name]
        
        feature = [feature for feature, codes in feature_defs.items()
                   if interaction[0] in codes]

        if 'hbond' in feature:
            feature.remove('hbond')
        assert len(feature) == 1
        feature = feature[0]

        if 'hbond_' in feature:
            feature = feature.replace('hbond_', '')
        return '{}{}:{}'.format(name, num, feature)
    
    def get_poses(self, cluster):
        return {l:self.ligands[l].poses[p] for l,p in cluster.items()}

    def get_rmsd(self, cluster):
        tmp = [self.ligands[l].poses[p].rmsd for l,p in cluster.items()]
        return np.mean([r for r in tmp if r is not None])

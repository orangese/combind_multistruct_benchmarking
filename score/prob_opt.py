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
    gc50 (float): Glide score for which 50% of poses are correct.

    ligands ([containers.Ligand,])
    ligand_names ([str, ])
    xtal (set([str,]))

    mcss (mcss.MCSSController)
    shape (mcss.ShapeController)

    single (np.array, # ligands x 2)
    pair (np.array, # ligands x # ligands x max_poses x max_poses)
    corr (np.array, # ligands x # ligands)
    """
    def __init__(self, stats, features, max_poses, alpha, gc50):
        self.stats = stats
        self.features = features
        self.max_poses = max_poses
        self.alpha = float(alpha)
        self.gc50 = float(gc50)

        self.ligands = []
        self.ligand_names = []
        self.xtal = set()

        self.mcss = None
        self.shape = None

        self.single = None
        self.pair = None
        self.corr = None

    def set_ligands(self, ligands, mcss, shape, xtal=set(), matrix=True):
        ligands = list(ligands.items())
        self.ligand_names = [x[0] for x in ligands]
        self.ligands = [x[1] for x in ligands]

        self.mcss = mcss
        self.shape = shape
        self.xtal = xtal

        if matrix:
            self.single = self._get_single()
            self.pair = self._get_pair()
            self.corr = self._get_corr()

    def _get_single(self):
        single = np.zeros((len(self.ligands), self.max_poses))-20
        for i, ligand in enumerate(self.ligands):
            for j, pose in enumerate(ligand.poses[:self.max_poses]):
                if self.gc50 == float('inf'):
                    x = -self.alpha*pose.gscore
                else:
                    x = -self.alpha*(pose.gscore-self.gc50)
                single[i, j] = x

                if ligand.ligand in self.xtal:
                    single[i, j] = 10
        return single

    def _get_pair(self):
        pair = np.zeros((len(self.ligands), len(self.ligands),
                         self.max_poses, self.max_poses))
        for i, ligand1 in enumerate(self.ligands):
            for j, ligand2 in enumerate(self.ligands[i+1:]):
                j += i+1
                lp = LigPair(ligand1, ligand2, self.features,
                             self.mcss, self.shape, self.max_poses)
                lp.init_pose_pairs()
                for r1, r2 in lp.pose_pairs:
                    lr = self._get_lr_pair(lp, r1, r2)
                    pair[i, j, r1, r2] = lr
                    pair[j, i, r2, r1] = lr
        return pair

    def _get_corr(self):
        corr = np.ones((len(self.ligands), len(self.ligands)))
        np.fill_diagonal(corr, 1)
        return corr

    def _get_lr_pair(self, lp, r1, r2):
        lr = 0
        for feature in self.features:
            x = lp.get_feature(feature, r1, r2)
            if x is None:
                continue
            if 'native' in self.stats:
                p_x_native  = self.stats['native'][feature](x)
                p_x         = self.stats['reference'][feature](x)
                lr += np.log(p_x_native) - np.log(p_x)
            else:
                lr += self.stats['feature'][0]*x + self.stats['feature'][1]
        return lr

    ###########################################################################
    def max_posterior(self, max_iterations, restart):
        """
        max_iterations (int): Maximum number of iterations to attempt before exiting.
        restart (int): Number of times to run the optimization
        """
        best_score, best_poses = -float('inf'), None
        for i in range(restart):
            if i == 0:
                poses = {lig: 0 for lig in self.ligand_names}
            else:
                poses = {lig: np.random.randint(self.max_poses)
                         for lig in self.ligand_names}

            poses = self.optimize_poses(poses, max_iterations)
            score = self.log_posterior(poses)
            if score > best_score:
                best_score = score
                best_poses = poses.copy()

            print(poses)
            print('run {}, score {}'.format(i, score))

        return best_poses

    def optimize_poses(self, poses, max_iterations):
        """
        poses ({ligand_name: current pose number, })
        max_iterations (int)
        """
        for _ in range(max_iterations):
            update = False
            for query in np.random.permutation(list(poses.keys())):
                iposes = self.poses_to_iposes(
                        {lig: pose for lig, pose in poses.items() if lig != query})
                iquery = self.ligand_names.index(query)
                best_pose = self.best_pose(iposes, iquery)
                if best_pose != poses[query]:
                    update = True
                    poses[query] = best_pose
            if not update:
                break
        return poses

    ###########################################################################

    def score_new_ligand(self, poses, ligand):
        best_pose, best_logp = -1, -float('inf')
        for pose in range(len(ligand.poses)):
            logp = self.partial_log_posterior(poses, ligand, pose)
            if logp >= best_logp:
                best_pose, best_logp = pose, logp
        return best_pose

    def partial_log_posterior(self, poses, query_ligand, query_pose):
        iposes = list(self.poses_to_iposes(poses).items())

        lr = -self.alpha*query_ligand.poses[query_pose].gscore

        for ligand, pose in iposes:
            lp = LigPair(self.ligands[ligand], query_ligand, self.features,
                         self.mcss, self.shape, self.max_poses)
            lr += self._get_lr_pair(lp, pose, query_pose) / len(poses)
        return lr

    ###########################################################################

    def best_pose(self, iposes, iquery):
        if self.gc50 == float('inf'):
            q = np.zeros((len(iposes), 2))
            q[:, 0] = 1.0
        else:
            q = self.get_prob(iposes)

        _corr = np.array([[self.corr[lig1, lig2] for lig2 in iposes] for lig1 in iposes])
        C = self.correlations(_corr, q).reshape(-1, 1)

        _pair = np.array([self.pair[iquery, lig, :, pose] for lig, pose in iposes.items()])

        p = np.vstack([self.single[iquery, :], np.zeros(self.single[iquery, :].shape)]).T
        p[:, 0] += np.sum(C*np.log(q[:, :1]*np.exp(_pair) + q[:, 1:]), axis=0)
        p = np.exp(p)
        p /= p.sum(axis=1, keepdims=True)
        return np.argmax(p[:, 0])

    def get_prob(self, iposes):
        single = np.array([self.single[l, p] for l, p in iposes.items()])
        pair = np.array([[self.pair[l1, l2, p1, p2] for l2, p2 in iposes.items()]
                          for l1, p1 in iposes.items()])
        corr = np.array([[self.corr[l1, l2] for l2 in iposes] for l1 in iposes])
        return self.message_passing(single, pair, corr)

    def message_passing(self, single, pair, corr, max_iter=100, tol=0.0001):
        q_old = np.vstack([single, np.zeros(single.shape)]).T
        q_old = np.exp(q_old)
        q_old /= q_old.sum(axis=1, keepdims=True)
        for _ in range(max_iter):
            q = np.vstack([single, np.zeros(single.shape)]).T
            C = self.correlations(corr, q_old)
            q[:, 0] += np.sum(C*np.log(q_old[:, 0]*np.exp(pair) + q_old[:, 1]), axis=1)
            q = np.exp(q)
            q /= q.sum(axis=1, keepdims=True)
            if np.max(np.abs(q - q_old)) < tol:
                break
            q_old = q
        else:
            print('Message passing did not converge.') 
        return q

    def correlations(self, corr, q):
        C = corr*q[:, 0].reshape(1, -1)
        np.fill_diagonal(C, 1)
        return 1/C.sum(axis=1)

    ###########################################################################
    def log_posterior(self, poses):
        """
        Returns the log posterior for pose cluster.

        This should only be used for debugging purposes, as when optimizing we
        only need to consider terms that involve the ligand being considered.
        """
        iposes = list(self.poses_to_iposes(poses).items())

        lr = 0
        for lig, pose in iposes:
            lr += (len(poses)-1)*self.single[lig, pose]

        for i, (lig1, pose1) in enumerate(iposes):
            for lig2, pose2 in iposes[i+1:]:
                lr += self.pair[lig1, lig2, pose1, pose2]
        return lr
 
    def poses_to_iposes(self, poses):
        return {self.ligand_names.index(lig): pose for lig, pose in poses.items()}

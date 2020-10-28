import numpy as np
import os

def pad(x, shape1, shape2=0, C=1000):
    if len(x.shape) == 1:
        y = np.zeros(shape1)+C
        y[:x.shape[0]] = x[:shape1]
    elif len(x.shape) == 2:
        y = np.zeros((shape1, shape2))
        y[:x.shape[0], :x.shape[1]] = x[:shape1, :shape2]
    else:
        assert False
    return y

class PosePrediction:
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
    def __init__(self, ligands, raw, stats, xtal,
                 features, max_poses, alpha, gc50):
        self.ligands = ligands
        self.raw = raw
        self.stats = stats
        self.xtal = xtal
        self.features = features
        self.max_poses = self._get_max_poses(max_poses)
        self.alpha = float(alpha)
        self.gc50 = float(gc50)

        self.single = self._get_single()
        self.pair = self._get_pair()
        self.corr = self._get_corr()

    def _get_max_poses(self, max_poses):
        actual = max(len(x) for x in self.raw['gscore'].values())
        return min(max_poses, actual)

    def _get_single(self):
        single = []
        for ligand in self.ligands:
            gscore = self.raw['gscore'][ligand]
            if ligand in self.xtal:
                gscore[:] = -20.0
            single += [gscore]

        single = [pad(x, self.max_poses) for x in single]
        single = np.vstack(single)

        if self.gc50 != float('inf'):
            single -= self.gc50
        single = -self.alpha*single
        return single

    def _get_pair(self):
        pair = np.zeros((len(self.ligands), len(self.ligands),
                         self.max_poses, self.max_poses))
        for i, ligand1 in enumerate(self.ligands):
            for j, ligand2 in enumerate(self.ligands[i+1:]):
                j += i+1
                for feature in self.features:
                    stats = self.stats[feature]
                    raw = self.raw[feature][(ligand1, ligand2)]
                    energy = np.log(stats['native'](raw)) - np.log(stats['reference'](raw))
                    energy = pad(energy, self.max_poses, self.max_poses)
                    pair[i, j] += energy
                    pair[j, i] += energy.T
        return pair

    def _get_corr(self):
        corr = np.ones((len(self.ligands), len(self.ligands)))
        np.fill_diagonal(corr, 1)
        return corr

    ###########################################################################
    def max_posterior(self, max_iterations, restart):
        """
        max_iterations (int): Maximum number of iterations to attempt before exiting.
        restart (int): Number of times to run the optimization
        """
        if len(self.ligands) == 1:
            return {self.ligands[0]: 0}
        best_score, best_poses = -float('inf'), None
        for i in range(restart):
            if i == 0:
                poses = {lig: 0 for lig in self.ligands}
            else:
                poses = {lig: np.random.randint(self.max_poses)
                         for lig in self.ligands}

            poses = self.optimize_poses(poses, max_iterations)
            score = self.log_posterior(poses)
            if score > best_score:
                best_score = score
                best_poses = poses.copy()

            print(poses)
            print('run {}, score {}'.format(i, score))
        q = self.get_prob(self.poses_to_iposes(poses))
        C = self.correlations(self.corr, q)
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
                iquery = self.ligands.index(query)
                best_pose = self.best_pose(iposes, iquery)
                if best_pose != poses[query]:
                    update = True
                    poses[query] = best_pose
            if not update:
                break
        return poses

    def anneal_poses(self, poses, max_iterations):
        """
        poses ({ligand_name: current pose number, })
        max_iterations (int)
        """
        for T in np.logspace(2, -2, 11):
            print(T)
            for i in range(1000):
                print(i, 'of', len(self.ligands)*self.max_poses)
                for query in np.random.permutation(list(poses.keys())):
                    iposes = self.poses_to_iposes(
                            {lig: pose for lig, pose in poses.items() if lig != query})
                    iquery = self.ligands.index(query)
                    probs = self.get_probs(iposes, iquery)[:, 0]
                    probs = np.exp(np.log(probs)/T)
                    probs /= probs.sum()
                    poses[query] = np.random.choice(range(len(probs)), p=probs)
        return self.optimize_poses(poses, max_iterations)

    def best_pose(self, iposes, iquery):
        p = self.get_probs(iposes, iquery)
        return np.argmax(p[:, 0])

    def get_probs(self, iposes, iquery):
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
        return p

    def get_poses_prob(self, poses):
        iposes = self.poses_to_iposes(poses)
        probs = self.get_prob(iposes)[:, 0]

        ligands = [(self.ligands.index(lig), lig) for lig in poses]
        ligands = sorted(ligands)
        ligands = [lig[1] for lig in ligands]
        return {lig: prob for lig, prob in zip(ligands, probs)}

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
        return {self.ligands.index(lig): pose for lig, pose in poses.items()}

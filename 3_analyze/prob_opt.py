import numpy as np

from pairs import LigPair

class LigSet:
    def __init__(self, ligs, num_poses, features, mcss, t=1):
        self.ligs = ligs
        self.lig_ordering = sorted(ligs.keys())

        self.mcss = mcss
        self.features = features

        self.all_poses = {l:lig.poses[:num_poses] for l, lig in ligs.items()}
        self.priors = {l:self.init_prior(l, t) for l in ligs}

        self.lig_pairs = self.init_lig_pairs() 

    def init_prior(self, l, t):
        probs = [np.exp(-p.gscore/float(t)) for p in self.all_poses[l]]
        return [pr/sum(probs) for pr in probs]

    def init_lig_pairs(self):
        pairs = {}
        for i1, l1 in enumerate(self.lig_ordering):
            for i2, l2 in enumerate(self.lig_ordering[i1+1:]):
                pairs[(l1,l2)] = LigPair(self.ligs[l1], self.ligs[l2],
                    self.all_poses[l1], self.all_poses[l2], self.features, self.mcss)
        return pairs

    def get_poses(self, cluster):
        return {l:self.all_poses[l][p] for l,p in cluster.items()}

    def get_rmsd(self, cluster):
        return np.mean([self.all_poses[l][p].rmsd for l,p in cluster.items()])

class PredictStructs:
    def __init__(self, ligset, evidence, features):
        self.ligset = ligset
        self.ev = evidence
        self.features = features

        self.log_posterior_cache = {}

    def posterior(self, pose_cluster, sample_l=None, verbose=False, en_landscape=False):
        # pose cluster = ligand_name: integer_pose_index
        log_prob = 0
        if sample_l is not None:
            log_prob += np.log(self.ligset.priors[sample_l][pose_cluster[sample_l]])
        else:
            for l, pnum in pose_cluster.items():
                log_prob += np.log(self.ligset.priors[l][pnum])
                if verbose: print l, self.ligset.priors[l][pnum]

        if len(self.features) > 0:
            for (l1,l2), lp in self.ligset.lig_pairs.items():
                if sample_l is not None and sample_l != l1 and sample_l != l2: 
                    continue

                p1,p2 = pose_cluster[l1], pose_cluster[l2]
                prior = self.ligset.priors[l1][p1]*self.ligset.priors[l2][p2]
                if (l1,l2,p1,p2) not in self.log_posterior_cache:
                    pp = lp.pose_pairs[(p1,p2)]
                    pos = 0
                    for k, k_def in self.features.items():
                        k_val = pp.get_feature(k, k_def)
                        if k_val is None: continue
                        p_x_native = self.ev.evaluate(k, k_val, 'native')
                        p_x_nnative = self.ev.evaluate(k, k_val, 'decoy')
                        p_x = p_x_native*prior + p_x_nnative*(1-prior)
                        if p_x_native == 0:
                            #print k, p_x_native, p_x_nnative, prior
                            p_x_native = 10**(-8)
                        pos += np.log(p_x_native*prior/p_x)
                        if verbose: print l1, l2, k, round(k_val,2), np.log(prior*p_x_native/p_x)
                    self.log_posterior_cache[(l1,l2,p1,p2)] = pos
                log_prob += self.log_posterior_cache[(l1,l2,p1,p2)] - np.log(prior)

        if en_landscape:
            return log_prob, self.ligset.get_rmsd(pose_cluster)
        return log_prob, None

    def max_posterior(self, verbose=False, sampling=3, en_landscape=False):
        initial_cluster = {l:0 for l in self.ligset.ligs}
        max_sc, best_cluster, all_clusters = self.optimize(initial_cluster,sampling=sampling, en_landscape=en_landscape)

        if verbose:
            print 'cluster -1, prob = {}'.format(max_sc)

        # run 10 more times starting randomly
        for i in range(15):
            rand_cluster = {}
            for l in self.ligset.ligs:
                rand_cluster[l] = np.random.randint(len(self.ligset.all_poses[l]))

            new_sc, new_cluster, clusters = self.optimize(rand_cluster,sampling=sampling, en_landscape=en_landscape)
            all_clusters.extend(clusters)

            if new_sc > max_sc:
                max_sc, best_cluster = new_sc, new_cluster
            if verbose:
                print 'cluster {}, prob = {}'.format(i, new_sc)

        return best_cluster, all_clusters

    def optimize(self, init_cluster, sampling=3, en_landscape=False):
        # maximize objective by randomly sampling ligands
        time_since_update = 0
        pose_cluster = init_cluster

        all_clusters = []

        while time_since_update < sampling*(len(self.ligset.ligs))**2:
            time_since_update += 1

            rand_lig = np.random.choice(self.ligset.lig_ordering)
            sample_l = rand_lig
            if en_landscape: sample_l = None

            max_sc, rmsd = self.posterior(pose_cluster, sample_l=sample_l, en_landscape=en_landscape)
            best_p = pose_cluster[rand_lig]
            if en_landscape: all_clusters.append((max_sc, rmsd))
            
            for p in range(len(self.ligset.all_poses[rand_lig])):
                pose_cluster[rand_lig] = p
                new_score, new_rmsd = self.posterior(pose_cluster, sample_l=sample_l, en_landscape=en_landscape)
                if new_score > max_sc:
                    max_sc, best_p = new_score, p
                    time_since_update = 0
                if en_landscape: all_clusters.append((new_score, new_rmsd))

            pose_cluster[rand_lig] = best_p

        return self.posterior(pose_cluster)[0], pose_cluster, all_clusters
            





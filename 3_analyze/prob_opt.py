import numpy as np

from pairs import LigPair

class LigSet:
    def __init__(self, docking_st, evidence, features, max_poses, t):
        self.docking_st = docking_st
        self.ev = evidence
        self.features = features
        self.max_poses = max_poses
        self.num_poses = {}
        self.t = float(t)

        self.log_prior_cache = {}
        self.log_posterior_cache = {}

        self.lig_pairs = {}

    def log_prior(self, l, r):
        if l not in self.log_prior_cache:
            t,n=self.t,self.n(l)
            probs = [np.exp(-p.gscore/t) for p in self.docking_st.ligands[l].poses[:n]]
            self.log_prior_cache[l] = [np.log(pr/sum(probs)) for pr in probs]
        return self.log_prior_cache[l][r]

    def get_lig_pair(self, l1, l2):
        if (l1,l2) in self.lig_pairs:
            return self.lig_pairs[(l1,l2)]#.pose_pairs[(p1,p2)]
        if (l2,l1) in self.lig_pairs:
            return self.lig_pairs[(l2,l1)]#.pose_pairs[(p2,p1)]

        self.lig_pairs[(l1,l2)] = LigPair(self.docking_st.ligands[l1], self.docking_st.ligands[l2],
                                          self.features, self.docking_st.mcss, self.max_poses, self.max_poses)

        return self.lig_pairs[(l1,l2)]#.pose_pairs[(p1,p2)]

    def n(self, l):
        if l not in self.num_poses:
            self.num_poses[l] = min(self.max_poses, len(self.docking_st.ligands[l].poses))
        return self.num_poses[l]

    def get_poses(self, cluster):
        return {l:self.docking_st.ligands[l].poses[p] for l,p in cluster.items()}

    def get_rmsd(self, cluster):
        return np.mean([self.docking_st.ligands[l].poses[p].rmsd for l,p in cluster.items()])

class PredictStructs:
    def __init__(self, docking_st, evidence, features, num_poses, t):
        self.ligset = LigSet(docking_st, evidence, features, num_poses, t)

    def x(self, pose_cluster, k, k_def, sample_l=None, lig_order=None):
        if lig_order is None:
            lig_order = pose_cluster.keys()

        x = np.zeros( (len(lig_order), len(lig_order)) )
        log_p_l_given_x = np.zeros( (len(lig_order), len(lig_order)) )
        for i, l1 in enumerate(lig_order):
            for j, l2 in list(enumerate(lig_order))[i+1:]:
                if sample_l is not None and sample_l != l1 and sample_l != l2:
                    continue
                #pp = self.ligset.get_pose_pair(l1, pose_cluster[l1], l2, pose_cluster[l2])
                lp = self.ligset.get_lig_pair(l1, l2)
                k_val = lp.get_feature(k, pose_cluster[l1],pose_cluster[l2])
                if k_val is None: continue
                assert k_val <= 1 and k_val >= 0, '{} {}'.format(fname, kval)
                x[i,j] = k_val
         
                prior = np.exp(self.ligset.log_prior(l1,pose_cluster[l1])+self.ligset.log_prior(l2,pose_cluster[l2]))
       
                p_x_native = self.ligset.ev.evaluate(k, k_val, 1)
                p_x_nnative = self.ligset.ev.evaluate(k, k_val, 0)

                # option A:
                p_x = p_x_native*prior + p_x_nnative*(1-prior)

                # option B:
                #p_x = self.ligset.ev.evaluate(k, k_val, -1)               

                log_p_l_given_x[i,j] = np.log(p_x_native/p_x)
                if p_x == 0:
                    print l1,l2,k, k_val
                    print p_x_native, p_x

        return x, log_p_l_given_x

    def joint_posterior(self, pose_cluster, sample_l=None, verbose=False, en_landscape=False):
        # pose cluster = ligand_name: integer_pose_index
        
        log_prob = sum([self.ligset.log_prior(l,p) for l,p in pose_cluster.items()])
        for k,k_def in self.ligset.features.items():
            x_k, log_p_l_given_x_k = self.x(pose_cluster, k, k_def, sample_l=sample_l)
            log_prob += np.sum(log_p_l_given_x_k)

        if en_landscape:
            return log_prob, self.ligset.get_rmsd(pose_cluster)
        return log_prob, None

    def max_posterior(self, ligands, verbose=False, sampling=3, en_landscape=False):
        initial_cluster = {l:0 for l in ligands}
        max_sc, best_cluster, all_clusters = self.optimize(initial_cluster,sampling=sampling, en_landscape=en_landscape)

        if verbose:
            print 'cluster -1, score {}'.format(max_sc)#, max_sc#, rmsd:', self.joint_posterior(best_cluster, en_landscape=True)

        # run 10 more times starting randomly
        for i in range(15):
            rand_cluster = {}
            for l in ligands:
                rand_cluster[l] = np.random.randint(self.ligset.num_poses[l])

            new_sc, new_cluster, clusters = self.optimize(rand_cluster,sampling=sampling, en_landscape=en_landscape)
            all_clusters.extend(clusters)

            if new_sc > max_sc:
                max_sc, best_cluster = new_sc, new_cluster
            if verbose:
                print 'cluster {}, score {}'.format(i,new_sc)# rmsd:'.format(i), self.joint_posterior(pose_cluster, en_landscape=True)

        return best_cluster, all_clusters

    def optimize(self, init_cluster, sampling=3, en_landscape=False):
        # maximize objective by randomly sampling ligands
        time_since_update = 0
        pose_cluster = init_cluster

        all_clusters = []

        while time_since_update < sampling*(len(init_cluster.keys()))**2:
            time_since_update += 1

            rand_lig = np.random.choice(init_cluster.keys())
            sample_l = rand_lig
            if en_landscape: sample_l = None

            max_sc, rmsd = self.joint_posterior(pose_cluster, sample_l=sample_l, en_landscape=en_landscape)
            best_p = pose_cluster[rand_lig]
            if en_landscape: all_clusters.append((max_sc, rmsd))
            
            for p in range(self.ligset.num_poses[rand_lig]):
                pose_cluster[rand_lig] = p
                new_score, new_rmsd = self.joint_posterior(pose_cluster, sample_l=sample_l, en_landscape=en_landscape)
                if new_score > max_sc:
                    max_sc, best_p = new_score, p
                    time_since_update = 0
                if en_landscape: all_clusters.append((new_score, new_rmsd))

            pose_cluster[rand_lig] = best_p

        return self.joint_posterior(pose_cluster)[0], pose_cluster, all_clusters
            
    def score_query(self, query, optimal_cluster):
        eval_cluster = {l:p for l,p in optimal_cluster.items()}
        scores = []
        for p in range(self.ligset.n(query)):
            eval_cluster[query] = p
            scores.append(self.joint_posterior(eval_cluster, sample_l=query)[0])
        return scores




import numpy as np

class ScoreQuery:
    def __init__(self, query_ligand, helper_ligands, all_lig, data_const, max_poses, w, cache={}):

        # two options:
        # 1 - pass in query ligand to score all poses of
        # 2 - don't pass in query to find top cluster among l_i

        self.l_q = query_ligand   # string
        self.l_i = helper_ligands # list of strings
        self.all_lig = all_lig    # { lig_id : Ligand}

        self.w = w
        self.d_const = data_const

        self.max_poses = max_poses

        self.cache = cache

    def score_all(self):
        self.num_poses = {l: min(self.max_poses, len(self.all_lig[l].poses)) for l in self.l_i}
        if self.l_q is not None:
            self.num_poses[self.l_q] = min(self.max_poses, len(self.all_lig[self.l_q].poses))
            self.pose_neighbors = {i:self.score(i) for i in range(self.num_poses[self.l_q])}
            self.pose_scores = {i: self.objective(cluster) for i, cluster in self.pose_neighbors.items()}
        else:
            self.top_cluster = self.score(None)

    def load(self, p_neighbors):
        self.pose_neighbors = p_neighbors
        self.pose_scores = {i: self.objective(cluster) for i, cluster in self.pose_neighbors.items()}

    def score(self, query_pose):
        # return map of ligands to poses that
        # maximize the objective

        # run optimization starting from glide poses
        initial_cluster = {}
        if query_pose is not None:
            initial_cluster = {self.l_q : self.all_lig[self.l_q].poses[query_pose]}
        for l in self.l_i:
            initial_cluster[l] = self.all_lig[l].poses[0]

        max_sc, best_cluster = self.optimize(initial_cluster)

        # run 10 more times starting randomly
        for i in range(10):
            rand_cluster = {}
            if query_pose is not None:
                rand_cluster = {self.l_q : self.all_lig[self.l_q].poses[query_pose]}
            for l in self.l_i:
                rand_cluster[l] = self.all_lig[l].poses[np.random.randint(self.num_poses[l])]
            new_sc, new_cluster = self.optimize(rand_cluster)
            if new_sc > max_sc:
                max_sc, best_cluster = new_sc, new_cluster

        return best_cluster

    def optimize(self, init_cluster, sampling=2):
        # maximize objective by randomly sampling ligands

        time_since_update = 0
        pose_cluster = init_cluster

        while time_since_update < sampling*(len(self.l_i))**2:
            time_since_update += 1

            rand_lig = np.random.choice(self.l_i)

            max_score, best_pose = self.objective(pose_cluster, l_sample=rand_lig), pose_cluster[rand_lig]
            for p in range(self.num_poses[rand_lig]):
                pose_cluster[rand_lig] = self.all_lig[rand_lig].poses[p]
                new_score = self.objective(pose_cluster, l_sample=rand_lig)
                if new_score > max_score:
                    max_score, best_pose = new_score, self.all_lig[rand_lig].poses[p]
                    time_since_update = 0
            pose_cluster[rand_lig] = best_pose

        return self.objective(pose_cluster), pose_cluster

    def objective(self, cluster, l_sample=None):
        sc = 0
        all_ligs = cluster.keys()
        for i1, l1 in enumerate(all_ligs):

            if l_sample is None:
                other_ligs = all_ligs[i1+1:]
            elif l1 != l_sample:
                continue
            else:
                other_ligs = [l for l in all_ligs if l != l1]

            sc += -1*cluster[l1].gscore

            for l2 in other_ligs:
                sc += self.d_const*self.score_pair(l1, cluster[l1], l2, cluster[l2])

        return sc/float(len(self.l_i)+1)

    def score_pair(self, l1, p1, l2, p2):
        if (l1, p1.rank, l2, p2.rank) in self.cache:
            return self.cache[(l1, p1.rank, l2, p2.rank)]
        elif (l2, p2.rank, l1, p1.rank) in self.cache:
            return self.cache[(l2, p2.rank, l1, p1.rank)]

        p1.weight_fp(self.w)
        p2.weight_fp(self.w)

        sc = 0
        #print p1.fp.keys()
        #print p2.fp.keys()
        for r in p1.fp:
            if r in p2.fp:
                sc += (p1.fp[r]*p2.fp[r])**0.5

        self.cache[(l1, p1.rank, l2, p2.rank)] = sc
        return sc 










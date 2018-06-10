import numpy as np

k_defs = {
    'mcss':[],
    'hbond':[2,3],
    'sb1':[0],
    'sb2':[1],
    'sb3':[4],
    #'pipi1':[5],
    'pipi':[6],
    'picat':[7,8],
    #'c1':[10],
    'contact':[11]
}

class LigPair:
    def __init__(self, l1, l2, features, mcss, max_poses, normalize_fp):
        self.l1 = l1
        self.l2 = l2
        self.max_poses = max_poses
        self.mcss = mcss
        self.k_list = features
        self.normalize_fp = normalize_fp

        self.pose_pairs = {}
        self.feat_map = None

        # Only need to consider all pose pairs if normalizing
        # Otherwise compute them lazily
        if self.normalize_fp:
            self.init_pose_pairs()
            self._init_feat_map()

    def get_feature(self, k, r1, r2):
        """
        Compute the feature value for poses r1 and r2 for feature 'f_name'.
        Normalize by the maximum feature value for any pose pair
        if self.normalize_fp.
        """
        pp = self.get_pose_pair(r1, r2)
        x_k = pp.get_feature(k)
        if x_k is None: return None # ignore feature
        
        # Normalization factor is 1 if unnormalized, otherwise
        # it is the maximum feature value for the pair
        mi, ma = 0.0, 1.0
        if self.normalize_fp:
            mi,ma = self.feat_map[k]
            if mi == ma: return None

        if k == 'mcss': return 1 - x_k/ma
        return x_k/ma

    def get_pose_pair(self, r1, r2):
        if (r1, r2) not in self.pose_pairs:
            l1, l2 = self.l1.lig_id, self.l2.lig_id
            mcss_score = None
            if (l1,l2) in self.mcss:
                mcss_score = self.mcss[(l1,l2)][(r1,r2)]
            if (l2,l1) in self.mcss:
                mcss_score = self.mcss[(l2,l1)][(r2,r1)]

            p1, p2 = self.l1.poses[r1], self.l2.poses[r2]
            assert p1.rank == r1, "Pose indices don't match rank"
            assert p2.rank == r2, "Pose indices don't match rank"
            self.pose_pairs[(r1,r2)] = PosePair(p1,p2,mcss_score)
        return self.pose_pairs[(r1, r2)]

    def init_pose_pairs(self):
        # r1 and r2 are ranks (glidescore indices)
        for r1 in range(min(len(self.l1.poses), self.max_poses)):
            for r2 in range(min(len(self.l2.poses), self.max_poses)):
                self.get_pose_pair(r1, r2)

    def _init_feat_map(self):
        self.feat_map = {k:(float('inf'),-float('inf')) for k in self.k_list}
        for key,pp in self.pose_pairs.items():
            for f,(minval,maxval) in self.feat_map.items():
                pp_x = pp.get_feature(f)
                self.feat_map[f] = (min(minval,pp_x),max(maxval,pp_x))


class PosePair:
    def __init__(self, p1, p2, mcss_score):
        self.p1 = p1
        self.p2 = p2
        self.x = {'mcss':mcss_score}

    def correct(self):
        if self.p1.rmsd <= 2 and self.p2.rmsd <= 2:
            return 1
        return 0

    def get_feature(self, k):
        if k not in self.x:
            self.x[k] = self.interaction_similarity(k)
        return self.x[k]

    def interaction_similarity(self, k):
        assert k in k_defs, k
        sc = 0 
        for i,r in self.p1.fp:
            if i in k_defs[k] and (i,r) in self.p2.fp:
                sc += (self.p1.fp[(i,r)]*self.p2.fp[(i,r)])**0.5
        return sc

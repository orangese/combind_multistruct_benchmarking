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

        self.pose_pairs = self.init_pose_pairs()
        self.feat_map = self._init_feat_map()

    def get_feature(self, k, r1, r2):
        x_k = self.pose_pairs[(r1,r2)].get_feature(k)

        mi,ma = self.feat_map[k]
        if mi == ma or x_k is None: return None

        if not self.normalize_fp: ma = 1
        if k == 'mcss' and self.normalize_fp: return 1 - x_k/ma
        return x_k/ma

    def init_pose_pairs(self):
        pairs = {}
        for r1 in range(min(len(self.l1.poses), self.max_poses)):
            for r2 in range(min(len(self.l2.poses), self.max_poses)):
                mcss_score = self.mcss.get_rmsd(self.l1.lig_id, self.l2.lig_id, r1,  r2)
                p1, p2 = self.l1.poses[r1], self.l2.poses[r2]
                assert p1.rank == r1, "Pose indices don't match rank"
                assert p2.rank == r2, "Pose indices don't match rank"
                pairs[(r1,r2)] = PosePair(p1,p2,mcss_score)
        return pairs

    def _init_feat_map(self):
        feat_map = {k:(float('inf'),-float('inf')) for k in self.k_list}
        for key,pp in self.pose_pairs.items():
            for f,(minval,maxval) in feat_map.items():
                pp_x = pp.get_feature(f)
                feat_map[f] = (min(minval,pp_x),max(maxval,pp_x))
        return feat_map

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

import numpy as np

class LigPair:
    def __init__(self, l1, l2, features, mcss, max_poses, normalize_fp):
        self.l1 = l1
        self.l2 = l2
        self.max_poses = max_poses
        self.mcss = mcss
        self.features = features
        self.normalize_fp = normalize_fp

        self.pose_pairs = {}
        self.feat_map = None

        # Only need to consider all pose pairs if normalizing
        # Otherwise compute them lazily
        if self.normalize_fp:
            self.init_pose_pairs()
            self._init_feat_map()

    def get_feature(self, f_name, r1, r2):
        """
        Compute the feature value for poses r1 and r2 for feature 'f_name'.
        Normalize by the maximum feature value for any pose pair
        if self.normalize_fp.
        """
        pp = self.get_pose_pair(r1, r2)
        pp_x = pp.get_feature(f_name,self.features[f_name])
        if pp_x is None: return None # ignore feature
        
        # Normalization factor is 1 if unnormalized, otherwise
        # it is the maximum feature value for the pair
        mi, ma = 0.0, 1.0
        if self.normalize_fp:
            mi,ma = self.feat_map[f_name]
            if mi == ma: return None

        if f_name == 'mcss': return 1 - pp_x/ma
        return pp_x/ma

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
        self.feat_map = {f:(float('inf'),-float('inf')) for f in self.features}
        for key,pp in self.pose_pairs.items():
            for f,(minval,maxval) in self.feat_map.items():
                pp_x = pp.get_feature(f,self.features[f])
                self.feat_map[f] = (min(minval,pp_x),max(maxval,pp_x))


class PosePair:
    def __init__(self, p1, p2, mcss_score):
        self.p1 = p1
        self.p2 = p2
    
        self.x = {'mcss':None}
        if mcss_score is not None:
            self.x['mcss'] = round(mcss_score,2)

    def correct(self):
        if self.p1.rmsd <= 2 and self.p2.rmsd <= 2:
            return 1
        return 0

    def get_feature(self, feat_name, feat_def):
        if feat_name in self.x:
            return self.x[feat_name]
        self.x[feat_name] = self.interaction_similarity(feat_def)
        return self.x[feat_name]

    def interaction_similarity(self, i_list):
        sc = 0 
        for i,r in self.p1.fp:
            if i in i_list and (i,r) in self.p2.fp:
                sc += (self.p1.fp[(i,r)]*self.p2.fp[(i,r)])**0.5
        return round(sc,2)

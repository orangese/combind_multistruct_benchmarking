import numpy as np

class LigPair:
    def __init__(self, l1, l2, features, mcss, max_r1=100, max_r2=100):
        self.l1 = l1
        self.l2 = l2

        self.max_r1 = max_r1
        self.max_r2 = max_r2

        self.mcss = mcss
        self.features = features

        self.pose_pairs = self.init_pose_pairs()
        self.feat_map = self.init_feat_map()

    def init_pose_pairs(self):
        # r1 and r2 are ranks (glidescore indices)
        pairs = {}
        for p1 in self.l1.poses[:self.max_r1]:
            for p2 in self.l2.poses[:self.max_r2]:#r2_list:
                l1, l2 = self.l1.lig_id, self.l2.lig_id
                #p1, p2 = self.l1.poses[r1], self.l2.poses[r2]
                mcss_score = None
                if (l1,l2) in self.mcss:
                    mcss_score = self.mcss[(l1,l2)][(p1.rank,p2.rank)]
                if (l2,l1) in self.mcss:
                    mcss_score = self.mcss[(l2,l1)][(p2.rank,p1.rank)]
                pairs[(p1.rank,p2.rank)] = PosePair(p1,p2,mcss_score)
        return pairs

    def init_feat_map(self):
        feat_map = {f:(float('inf'),-float('inf')) for f in self.features}
        for key,pp in self.pose_pairs.items():
            for f,(minval,maxval) in feat_map.items():
                pp_x = pp.get_feature(f,self.features[f])
                feat_map[f] = (min(minval,pp_x),max(maxval,pp_x))
        return feat_map

    def get_feature(self, f_name, r1, r2):
        mi,ma = self.feat_map[f_name]
        pp_x = self.pose_pairs[(r1,r2)].get_feature(f_name,self.features[f_name])
        if mi == ma or pp_x is None: return None # ignore feature
        if f_name == 'mcss': 
            return 1 - pp_x/ma
        return pp_x/ma #(pp_x - mi)/(ma - mi)

class PosePair:
    def __init__(self, p1, p2, mcss_score): #, active_features):
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




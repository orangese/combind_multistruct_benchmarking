import numpy as np

class LigPair:
    def __init__(self, l1, l2, features, mcss):
        self.l1 = l1
        self.l2 = l2

        self.mcss = mcss
        self.active_features = [f for f in features if self.active(f)]

        self.pose_pairs = {}
        
    def active(self, f):
        if f == 'mcss':
            l1,l2 = self.l1.lig_id, self.l2.lig_id
            if (l1,l2) in self.mcss: return True
            if (l2,l1) in self.mcss: return True
            return False
        if f in ['hbond','contact']: return True
        if f == 'pipi':
            return self.l1.ring and self.l2.ring
        if f == 'sb':
            return self.l1.chrg and self.l2.chrg

    def get_pose_pair(self, r1, r2):
        # r1 and r2 are ranks (glidescore indices)
        if (r1,r2) not in self.pose_pairs:
            l1, l2 = self.l1.lig_id, self.l2.lig_id
            p1, p2 = self.l1.poses[r1], self.l2.poses[r2]
            mcss_score = None
            if (l1,l2) in self.mcss:
                mcss_score = self.mcss[(l1,l2)][(r1,r2)]
            if (l2,l1) in self.mcss:
                mcss_score = self.mcss[(l2,l1)][(r2,r1)]
            self.pose_pairs[(r1,r2)] = PosePair(p1,p2,mcss_score,self.active_features)
        return self.pose_pairs[(r1,r2)]

    def all_pose_pairs(self, r1_list, r2_list):
        return [self.get_pose_pair(r1,r2) for r1 in r1_list for r2 in r2_list] 


class PosePair:
    def __init__(self, p1, p2, mcss_score, active_features):
        self.p1 = p1
        self.p2 = p2

        self.active_features = set(active_features)
    
        self.x = {}
        if mcss_score is not None:
            self.x['mcss'] = round(mcss_score,2)
    
    def native(self):
        return self.p1.rmsd <= 2 and self.p2.rmsd <= 2

    def correct(self):
        if self.p1.rmsd <= 2 and self.p2.rmsd <= 2:
            return 1
        return 0

    def get_feature(self, feat_name, feat_def):
        if feat_name in self.x:
            return self.x[feat_name]
        if feat_name not in self.active_features: return None
        self.x[feat_name] = self.interaction_similarity(feat_def)
        return self.x[feat_name]

    def interaction_similarity(self, i_list):
        sc = 0 
        for i,r in self.p1.fp:
            if i in i_list and (i,r) in self.p2.fp:
                sc += (self.p1.fp[(i,r)]*self.p2.fp[(i,r)])**0.5
        return round(sc,2)




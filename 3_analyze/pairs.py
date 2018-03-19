import numpy as np

class LigPair:
    def __init__(self, l1, l2, p1_list, p2_list, features, mcss):
        self.l1 = l1
        self.l2 = l2
        self.p1_list = p1_list
        self.p2_list = p2_list

        self.mcss = mcss
        self.active_features = [f for f in features if self.active(f)]

        self.pose_pairs = self.initialize_pp()
        
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

    def initialize_pp(self):
        l1, l2 = self.l1.lig_id, self.l2.lig_id
        tr = {}
        for p1 in self.p1_list:
            for p2 in self.p2_list:
                mcss_score = None
                if (l1,l2) in self.mcss:
                    mcss_score = self.mcss[(l1,l2)][(p1.rank,p2.rank)]
                if (l2,l1) in self.mcss:
                    mcss_score = self.mcss[(l2,l1)][(p2.rank,p1.rank)]
                tr[(p1.rank,p2.rank)] = PosePair(p1,p2,mcss_score,self.active_features)
        return tr

    def all_pose_pairs(self):
        return [pp for ij,pp in self.pose_pairs.items()]    


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




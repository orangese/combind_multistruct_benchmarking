import numpy as np

class LigPair:
    """
    Computes and stores overlap scores for a ligand pair.
    """
    def __init__(self, l1, l2, features, mcss, max_poses,
                 mode='maxoverlap', native_thresh=2.0):
        self.l1 = l1
        self.l2 = l2
        self.max_poses = max_poses
        self.mcss = mcss
        self.features = features
        self.mode = mode
        self.native_thresh = native_thresh

        self.pose_pairs = self._init_pose_pairs()
        self.feat_map = self._init_feat_map()

    def get_feature(self, feature, rank1, rank2):
        """
        Returns the overlap score for feature for the pose pair
        consisting of poses rank1 and rank2.
        """
        if feature == 'mcss':
            return self.pose_pairs[(rank1,rank2)].mcss_score

        if 'exclude' in self.mode or 'maxoverlap' in self.mode:
            minimum, maximum = self.feat_map[feature]
            if maximum == 0:
                return None

        if self.mode == 'maxoverlap' or (self.mode == 'hybrid' and feature == 'contact'):
            feature_value = self.pose_pairs[(rank1,rank2)].overlap(feature)
            minimum, maximum = self.feat_map[feature]
            return feature_value / max(maximum, 1.0)
        elif 'tanimoto-maxoverlap' in self.mode:
            minimum, maximum = self.feat_map[feature]
            return self.pose_pairs[(rank1, rank2)].tanimoto(feature, maxoverlap=maximum)
        else:
            return self.pose_pairs[(rank1, rank2)].tanimoto(feature)

    def _init_pose_pairs(self):
        """
        Create pose pairs for up to the top self.max_poses poses.
        """
        pairs = {}
        for rank1 in range(min(len(self.l1.poses), self.max_poses)):
            for rank2 in range(min(len(self.l2.poses), self.max_poses)):
                if 'mcss' in self.features:
                    mcss_score = self.mcss.get_rmsd(self.l1.ligand, self.l2.ligand,
                                                    rank1, rank2)
                else:
                    mcss_score = None
                pairs[(rank1, rank2)] = PosePair(self.l1.poses[rank1],
                                                 self.l2.poses[rank2],
                                                 mcss_score,
                                                 self.features,
                                                 self.native_thresh)
        return pairs

    def _init_feat_map(self):
        """
        Compute the maximum and minimum value of each feature.
        """
        feat_map = {feature: (float('inf'), -float('inf'))
                    for feature in self.features}
        for pp in self.pose_pairs.values():
            for feature in self.features:
                if feature == 'mcss':
                    pp_x = pp.mcss_score
                else:
                    pp_x = pp.overlap(feature)
                if pp_x is None: continue
                feat_map[feature] = (min(feat_map[feature][0], pp_x),
                                     max(feat_map[feature][1], pp_x))
        return feat_map

class PosePair:
    """
    Computes overlap scores for a pair of poses.
    """
    def __init__(self, pose1, pose2, mcss_score, features, native_thresh):
        self.pose1 = pose1
        self.pose2 = pose2
        self.mcss_score = mcss_score
        self.features = features
        self.native_thresh = native_thresh

    def correct(self):
        return int(max(self.pose1.rmsd,self.pose2.rmsd) <= self.native_thresh)

    def overlap(self, feature):
        overlap = 0
        for (i, r) in self.pose1.fp:
            if i in self.features[feature] and (i, r) in self.pose2.fp:
                overlap += self._residue_level_overlap(self.pose1.fp[(i,r)],
                                                       self.pose2.fp[(i,r)])
        return overlap

    def tanimoto(self, feature, pseudo_hits=1, pseudo_misses=1, maxoverlap=float('inf')):
        overlap = pseudo_hits + self.overlap(feature)
        total = (2*pseudo_hits+pseudo_misses) + min(self._total(feature), 2*maxoverlap)
        return overlap / (total - overlap)

    def _residue_level_overlap(self, fp1, fp2):
        return (fp1*fp2)**0.5

    def _total(self, feature):
        total = 0
        for (i, r) in self.pose1.fp:
            if i in self.features[feature]:
                total += self.pose1.fp[(i,r)]

        for (i, r) in self.pose2.fp:
            if i in self.features[feature]:
                total += self.pose2.fp[(i,r)]
        return total

import numpy as np
from shared_paths import shared_paths, feature_defs

class LigPair:
    """
    Computes and stores overlap scores for a ligand pair.

    Importantly, this class normalizes the overlap scores
    by dividing by the maximum value.
    """
    def __init__(self, l1, l2, features, mcss, max_poses):
        self.l1 = l1
        self.l2 = l2
        self.max_poses = max_poses
        self.mcss = mcss
        self.features = features

        self.pose_pairs = self._init_pose_pairs()
        self.feat_map = self._init_feat_map()

    def get_feature(self, feature, rank1, rank2):
        """
        Returns the overlap score for feature for the pose pair
        consisting of poses rank1 and rank2.
        """
        feature_value = self.pose_pairs[(rank1,rank2)].get_feature(feature)
        minimum, maximum = self.feat_map[feature]
        
        if not maximum or feature_value is None:
            return None
        
        if feature == 'mcss':
            return feature_value
        return feature_value / max(maximum, 1.0)

    def get_gscores(self, rank1, rank2):
        pp = self.pose_pairs[(rank1,rank2)]
        return pp.pose1.gscore, pp.pose2.gscore

    def get_rmsds(self, rank1, rank2):
        pp = self.pose_pairs[(rank1,rank2)]
        return pp.pose1.rmsd, pp.pose2.rmsd

    def _init_pose_pairs(self):
        """
        Create pose pairs for up to the top self.max_poses poses.
        """
        pairs = {}
        for rank1 in range(min(len(self.l1.poses), self.max_poses)):
            for rank2 in range(min(len(self.l2.poses), self.max_poses)):
                if self.mcss:
                    mcss_score = self.mcss.get_rmsd(self.l1.ligand, self.l2.ligand, rank1, rank2)
                else:
                    mcss_score = None
                pose1, pose2 = self.l1.poses[rank1], self.l2.poses[rank2]
                pairs[(rank1,rank2)] = PosePair(pose1, pose2, mcss_score)
        return pairs

    def _init_feat_map(self):
        """
        Compute the maximum and minimum value of each feature.
        """
        feat_map = {feature: (float('inf'), -float('inf')) for feature in self.features}
        for key, pose_pair in self.pose_pairs.items():
            for f, (minval, maxval) in feat_map.items():
                pp_x = pose_pair.get_feature(f)
                if pp_x is None: continue
                feat_map[f] = (min(minval, pp_x), max(maxval, pp_x))
        return feat_map

class PosePair:
    """
    Computes and stores overlap scores for pose1 and pose2.

    This class defines the feature overlap scores as the sum of
    the geometric mean of the fingerprint values for each residue.
    """
    def __init__(self, pose1, pose2, mcss_score):
        self.pose1 = pose1
        self.pose2 = pose2
        self.features = {'mcss': mcss_score}

    def correct(self):
        """
        Returns 1 if both poses are at most 2 A RMSD from their
        crystallographic pose.
        """
        return int(    self.pose1.rmsd <= shared_paths['stats']['native_thresh']
                   and self.pose2.rmsd <= shared_paths['stats']['native_thresh'])

    def get_feature(self, feature):
        """
        Get the value of feature for this pose pair.
        """
        if feature not in self.features:
            self.features[feature] = self.interaction_similarity(feature)
        return self.features[feature]

    def interaction_similarity(self, feature):
        """
        Compute the overlap score for feature. Specifically,
        computes residue level scores by taking the geometric
        mean of the fingerprint values, then sums them to get
        a target level score.
        """
        assert feature in feature_defs, feature
        score = 0
        # (feature_index, residue)
        for (i, r) in self.pose1.fp:
            if i in feature_defs[feature] and (i, r) in self.pose2.fp:
                # Geometric mean
                score += (self.pose1.fp[(i,r)]*self.pose2.fp[(i,r)])**0.5
        return score

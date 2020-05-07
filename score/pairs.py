import numpy as np

class LigPair:
    """
    Computes and stores overlap scores for a ligand pair.
    """
    def __init__(self, l1, l2, features, mcss, max_poses, native_thresh=2.0):
        self.l1 = l1
        self.l2 = l2
        self.max_poses = max_poses
        self.mcss = mcss
        self.features = features
        self.native_thresh = native_thresh

        self.pose_pairs = {}

    def get_feature(self, feature, rank1, rank2):
        """
        Returns the overlap score for feature for the pose pair
        consisting of poses rank1 and rank2.
        """
        if feature not in self.features:
            return None

        key = (rank1, rank2)
        if key not in self.pose_pairs:
            self._init_pose_pair(rank1, rank2)

        if feature == 'mcss':
            return self.pose_pairs[key].mcss_score
        return self.pose_pairs[key].tanimoto(feature)

    def correct(self, rank1, rank2):
        key = (rank1, rank2)
        if key not in self.pose_pairs:
            self._init_pose_pair(rank1, rank2)

        return self.pose_pairs[key].correct()

    def _init_pose_pair(self, rank1, rank2):
        pose1 = self.l1.poses[rank1]
        pose2 = self.l2.poses[rank2]
        if 'mcss' in self.features:
            mcss_score = self.mcss.get_rmsd(self.l1.ligand, self.l2.ligand,
                                            pose1.rank, pose2.rank)
        else:
            mcss_score = None

        self.pose_pairs[(rank1, rank2)] = PosePair(pose1, pose2,
                                                   mcss_score,
                                                   self.features,
                                                   self.native_thresh)

    def init_pose_pairs(self):
        for rank1 in range(min(len(self.l1.poses), self.max_poses)):
            for rank2 in range(min(len(self.l2.poses), self.max_poses)):
                self._init_pose_pair(rank1, rank2)

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

    def tanimoto(self, feature, pseudo_hits=1, pseudo_misses=1):
        overlap = pseudo_hits + self.overlap(feature)
        total = (2*pseudo_hits+pseudo_misses) + self._total(feature)
        return overlap / (total - overlap)

    def overlap(self, feature):
        overlap = 0
        for (i, r) in self.pose1.fp:
            if i in self.features[feature] and (i, r) in self.pose2.fp:
                overlap += self._residue_level_overlap(self.pose1.fp[(i,r)],
                                                       self.pose2.fp[(i,r)])
        return overlap

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

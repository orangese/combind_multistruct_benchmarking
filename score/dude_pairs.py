import numpy as np
from shared_paths import shared_paths, feature_defs

class DUDEPDBLigPair:
    """
    Computes and stores overlap scores for a ligand pair between a DUDE ligand and a PDB ligand.

    Importantly, this class normalizes the overlap scores
    by dividing by the maximum value.

    Inputs:
    * dude_ligand (Ligand object)
    * pdb_ligand (Ligand object)
    * max_poses (int)
    * mcss (bool)
    * features (i.e. [feature (str) ...]) (list)
    """
    def __init__(self, dude_ligand, pdb_ligand, features, mcss, max_poses):
        self.dude_ligand = dude_ligand
        self.pdb_ligand = pdb_ligand
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
        ''' Returns glide scores of rank1 pose of the dude ligand and rank2 pose of the pdb ligand '''
        pp = self.pose_pairs[(rank1,rank2)]
        return pp.dude_pose.gscore, pp.pdb_pose.gscore

    def get_rmsds(self, rank1, rank2):
        pp = self.pose_pairs[(rank1,rank2)]
        return pp.dude_pose.rmsd, pp.pdb_pose.rmsd

    def _init_pose_pairs(self):
        """ Create PosePair's for the top pose of the DUDe ligand and up to the top self.max_poses 
        poses of the PDB ligand

        Returns:
        * pairs (dict): maps from (dude_rank (int), pdb_rank (int)) tuples, composed of the ranks of
        the two ligands -> PosePair object made from the two poses
        """
        pairs = {}
        # for rank1 in range(min(len(self.dude_ligand.poses), self.max_poses)):
        for rank1 in range(min(len(self.dude_ligand.poses), 1)):
            for rank2 in range(min(len(self.pdb_ligand.poses), self.max_poses)):
                if self.mcss:
                    mcss_score = self.mcss.get_rmsd(self.dude_ligand.ligand, self.pdb_ligand.ligand,
                                                    rank1, rank2)
                else:
                    mcss_score = None
                pairs[(rank1, rank2)] = DUDEPDBPosePair(self.dude_ligand.poses[rank1],
                                                 self.pdb_ligand.poses[rank2],
                                                 mcss_score)
        return pairs

    def _init_feat_map(self):
        """ Compute the maximum and minimum value of each feature.

        Returns:
        * feat_map (dict): maps feature (i.e. sb, mcss, etc.) -> (min_score (float), max_score (float)),
        where the two floats are the lowest and highest scores, respectively, observed of that feature type
        across all pose pairs of the two ligands
        """
        feat_map = {feature: (float('inf'), -float('inf')) for feature in self.features}
        for key, pose_pair in self.pose_pairs.items():
            for feature in self.features:
                pp_x = pose_pair.get_feature(feature)
                if pp_x is None: continue
                feat_map[feature] = (min(feat_map[feature][0], pp_x),
                                     max(feat_map[feature][1], pp_x))
        return feat_map

class DUDEPDBPosePair:
    """
    Computes and stores overlap scores for dude_pose and pdb_pose.

    This class defines the feature overlap scores as the sum of
    the geometric mean of the fingerprint values for each residue.
    """
    def __init__(self, dude_pose, pdb_pose, mcss_score):
        self.dude_pose = dude_pose
        self.pdb_pose = pdb_pose
        self.features = {'mcss': mcss_score}

    def is_DUDE_native(self):
        """
        Returns 1 if both poses are at most 2 A RMSD from their
        crystallographic pose.
        """
        return int(self.pdb_pose.rmsd <= shared_paths['stats']['native_thresh'])

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
        for (i, r) in self.dude_pose.fp:
            if i in feature_defs[feature] and (i, r) in self.pdb_pose.fp:
                # Geometric mean
                score += (self.dude_pose.fp[(i,r)]*self.pdb_pose.fp[(i,r)])**0.5
        return score

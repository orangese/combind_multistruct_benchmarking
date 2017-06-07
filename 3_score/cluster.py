class Cluster:
    def __init__(self, poses): #poses is a list of Pose objects
        self.floatsPerResidue = 6 #Could change if we change the fingerprint
        self.poses = poses
        self.residues = list(set().union(*[pose.fp.feats.keys() for pose in self.poses])) #Cache all of the residues
        self.numpyFPs = []
        self.score = None
        
        #Setup the numpy vectorized FPs
        for pose in self.poses:
            lengthFPVector = self.floatsPerResidue * len(self.residues)
            poseFP = np.zeros(lengthFPVector)
            rawPoseFP = pose.fp.feats
            
            for residueIndex, residue in enumerate(self.residues):
                if residue in rawPoseFP:
                    poseFP[residueIndex * self.floatsPerResidue:(residueIndex+1) * self.floatsPerResidue] = rawPoseFP[residue]
                    
            self.numpyFPs.append(poseFP)

    def get_individual_overlap(self, a):
        if np.sign(a[0]) == np.sign(a[1]): #If they are the same sign
            return np.minimum(np.abs(a[0]), np.abs(a[1]))
        return 0
    
    def overlap(self):
        totalOverlap = 0
        
        for indOne, fpOne in enumerate(self.numpyFPs):
            for indTwo, fpTwo in enumerate(self.numpyFPs):
                if indOne != indTwo:
                    newFP = np.vstack((fpOne, fpTwo))
                    totalOverlap += np.sum(np.apply_along_axis(self.get_individual_overlap, 0, newFP))
                
        return totalOverlap
    
    def overlap_per_pair(self):
        totalOverlap_ind = []
        ind_track = []
        
        for indOne, fpOne in enumerate(self.numpyFPs):
            for indTwo, fpTwo in enumerate(self.numpyFPs):
                if indOne != indTwo:
                    newFP = np.vstack((fpOne, fpTwo))
                    totalOverlap_ind.append(np.sum(np.apply_along_axis(self.get_individual_overlap, 0, newFP)))
                    ind_track.append((indOne, indTwo))
                
        return totalOverlap_ind, ind_track

    def __str__(self):
        return str(self.overlap())

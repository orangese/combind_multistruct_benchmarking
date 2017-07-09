import numpy as np

class FuzzyFingerPrint:
    def __init__(self, feats, pdb = None):
        self.pdb = pdb
        self.feats = feats
    
    @classmethod
    def compact_parser(cls, line, pdb=None, w=[10,10,10,1]):
        feats = {i[0] : map(float, i[1:]) for i in map(lambda x: x.split(','), line.split(':'))}

        w = np.array(w)
        for res in feats.keys():
            feats[res] = np.array([i**0.5 if i > 0 else -((-i)**0.5) for i in np.multiply(w,np.array(feats[res]))])

        return cls(feats,pdb)

import numpy as np

class FuzzyFingerPrint:
    def __init__(self, feats, pdb = None):
        self.pdb = pdb
        self.feats = feats
    
    @classmethod
    def compact_parser(cls, line, pdb=None, w=[10,10,10,0,1]):
        feats = {i[0] : map(float, i[1:]) for i in map(lambda x: x.split(','), line.split(':'))}

        new_feats = {}

        for r in feats:
            if int(r) > 1000 and str(int(r) - 1000) in feats:
                print r
                raise Exception
            elif int(r) > 1000:
                new_feats[str(int(r) - 1000)] = feats[r]
            else:
                new_feats[r] = feats[r]

        w = np.array(w)
        for res in new_feats.keys():
            new_feats[res] = np.array([i**0.5 if i > 0 else (-i)**0.5 for i in np.multiply(w,np.array(new_feats[res]))])

        return cls(new_feats,pdb)

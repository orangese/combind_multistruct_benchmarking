import numpy as np

class FuzzyFingerPrint:
    def __init__(self, feats, pdb = None):
        self.pdb = pdb
        self.feats = feats
    
    @classmethod
    def compact_parser(cls, line, pdb=None):
        feats = {i[0] : map(float, i[1:]) for i in map(lambda x: x.split(','), line.split(':'))}

        return cls(feats,pdb)

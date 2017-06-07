from GPCRdb_numbers.GPCRdb_numbers_aaNum import GPCRdbNumbers
import numpy as np

class PdbToUniprot:
    """
    Maps from PDBs to Uniprot IDs

    Only applies to class A GPCRs.
    If a PDB has multiple proteins, only the uniprot for the class A GPCR is recorded
    """
    def __init__(self, gpcrs, _pdbtosp = '/share/PI/rondror/docking/fingerprint/GPCRdb_numbers/pdbtosp.txt'):
        with open(_pdbtosp) as pdbtosp:
            self.uni = {}
            line = pdbtosp.readline().split()
            # Kinda sketchy
            while not line or line[0] != '101M': line = pdbtosp.readline().split()
            while line:
                    pdb, experiment = line[:2]
                    if experiment == 'X-ray':
                        for uniprot in [entry for entry in line[4:] if '_' in entry]:
                            if gpcrs.contains(uniprot):
                                self.uni[pdb] = uniprot
                    line = pdbtosp.readline().split()

    def get_sp(self, pdb):
        try:
            return self.uni[pdb.upper()]
        except:
            return 'not found'

    def contains(self, pdb):
        return pdb.upper() in self.uni

class FuzzyFingerPrint:
    def __init__(self, feats, pdb = None):
        self.pdb = pdb
        self.feats = feats
    
    @classmethod
    def compact_parser(cls, line, pdb=None):
        feats = {i[0] : map(float, i[1:]) for i in map(lambda x: x.split(','), line.split(':'))}

        return cls(feats,pdb)

     

class FingerPrint:
    INTERACTIONS = ('T-STACK', 'PI-PI', 'CAT-PI', 'SALT', 'HYDRO', 'DONOR', 'ACCEPT')
    GENERIC_MAP = GPCRdbNumbers()
    PDB_MAP = PdbToUniprot(GENERIC_MAP)

    def __init__(self, feats, pdb = None):
        self.pdb = pdb
        self.feats = feats

    @classmethod
    def _map_to_generic(cls, feats, pdb):
        uniprot = cls.PDB_MAP.get_sp(pdb)
        gen = {}
        for feat in feats:
            key, residue = cls.GENERIC_MAP.get_generic_num(uniprot, feat)
            if key == 'None': key = feat
            gen[key] = feats[feat]
        return gen

    def to_common(self, resi):
        if 'x' not in resi: return resi
        uniprot = self.PDB_MAP.get_sp(self.pdb)
        key, residue = self.GENERIC_MAP.get_aa_num(uniprot, resi)
        return key
        
    @classmethod
    def compact_parser(cls, line, pdb = None, map_to_generic = False):
        """
        Populate feats from compact ifp representation
        RESI,<CSV of interactions:RESI,<CSV of interactions>...
        Do not change residue ids
        """
        feats = {i[0] : map(float, i[1:]) for i in map(lambda x: x.split(','), line.split(':'))}
        if map_to_generic: feats = cls._map_to_generic(feats, pdb)
        return cls(feats, pdb)

    @classmethod
    def round(cls, string):
        val = int(float(string)*1)
        if val < -2: val = -3
        if val > 2: val = 3
        return val

    @classmethod
    def rounded_parser(cls, line, pdb = None, map_to_generic = False):
        """
        Populate feats from compact ifp representation
        RESI,<CSV of interactions:RESI,<CSV of interactions>...
        Do not change residue ids
        """
        feats = {i[0] : map(cls.round, i[1:]) for i in map(lambda x: x.split(','), line.split(':'))}
        if map_to_generic: feats = cls._map_to_generic(feats, pdb)
        return cls(feats, pdb)

    @classmethod
    def binana_parser(cls, path, pdb = None):
        """
        Populate self.feats from binana output.
        Map residue IDs to generic if applicable.
        """
        binana = BinanaParser(path)  #Creates a Binana_parser object
        feats = {}

        ##Getting the T-stacking features. Feature 0
        if len(binana.T_stacking_dict)!=0:
            for resi in binana.T_stacking_dict['receptor_resi']:
                location=resi[1]  #location of the residue
                if location not in feats: feats[location] = [0] * 7
                feats[location][0]=1

        #GETTING PI_PI stacking features. Feature 1
        if(len(binana.PiPi_stacking_dict))!=0:
            for resi in binana.PiPi_stacking_dict['receptor_resi']:
                location=resi[1]
                if location not in feats: feats[location] = [0] * 7
                feats[location][1]=1
            
            
        ##GETTING CATION-PI interactions. FEATURE 2  
        if(len(binana.CationPi_dict))!=0:
            for resi in binana.CationPi_dict['receptor_resi']:
                location=resi[1]
                if location not in feats: feats[location] = [0] * 7    
                feats[location][2]=1
            
        ##GETTING SALT-BRIDGES .         FEATURE 3
        if(len(binana.SaltBridge_dict)) !=0:
            for resi in binana.SaltBridge_dict['receptor_resi']:
                location=resi[1]
                if location not in feats: feats[location] = [0] * 7                
                feats[location][3]=1

        ##GETTING hydrophobic interactions .   FEATURE 4
        if(len(binana.hydrophobic_dict)) !=0:
            for resi in binana.hydrophobic_dict['receptor']:
                location=resi[1]
                if location not in feats: feats[location] = [0] * 7                
                feats[location][4]=1

        ##RECEPTOR is a H-bond donor.      FEATURE 5
        if (len(binana.H_bond_list)) !=0:
            for bond in binana.H_bond_list:
                if bond[0]=='RECEPTOR':   #ie the receptor is the donor
                    location=bond[1].strip('()').split('(')[1]
                    if location not in feats: feats[location] = [0] * 7                    
                    feats[location][5]=1

        ##RECEPTOR is a H-bond acceptor   FEATURE 6
        if(len(binana.H_bond_list)) !=0:
            for bond in binana.H_bond_list:
                if bond[0]=='LIGAND':
                    location=bond[1].strip('()').split('(')[1]
                    if location not in feats: feats[location] = [0] * 7    
                    feats[location][6]=1
        feats = cls._map_to_generic(feats, pdb) if cls.PDB_MAP.contains(pdb) else feats
        return cls(feats, pdb)

    def residues(self):
        return self.feats.keys()

    def vectorize(self, residues):
        return np.array([bit for resi in residues for bit in self.entry(resi)])

    def entry(self, residue):
        """
        Given a residue return FP for that position,
        If not present return zeros + 1, to indicate such
        """
        if residue in self.feats:
            return self.feats[residue] + [0]
        else:
            return [0] * 7 + [1]

    def overlap(self, other):
        return sum(sum(i*j for i, j in zip(self.entry(resi), other.entry(resi))) for resi in self.residues())

    def total(self):
        return sum(i for resi in self.residues() for i in self.entry(resi))

    def hydro_overlap(self, other):
        overlap = 0
        for resi in self.residues():
            hydro1 = any(self.entry(resi)[i] for i in self.HYDRO)
            hydro2 = any(other.entry(resi)[i] for i in self.HYDRO)
            overlap += hydro1*hydro2
            for j in range(len(self.entry(resi))):
                if j in self.HYDRO: continue
                overlap += self.entry(resi)[j] * other.entry(resi)[j]
        return overlap

    def hydro_footprint(self):
        return self.hydro_overlap(self)

    def tanimoto_coef(self, other):
        """
        http://www.sciencedirect.com/science/article/pii/S1471489216300674

        Comparison of binary fingerprints is most often performed using the Tanimoto
        coefficient (Tc) which is the number of common bits in the two fingerprints 
        divided by the number of bits present in at least one of the fingerprints. 
        Tc ranges from 0 for dissimilar binding interactions to 1 for identical interactions [13]. 
        0.6 is generally considered as a minimum cutoff for binding mode similarity [15] 
        """
        overlap = self.overlap(other)
        flatten = lambda x: [item for sublist in x for item in sublist]
        total_bits = sum(flatten(self.feats.values()) + flatten(other.feats.values()))
        union_size = total_bits - overlap
        return overlap / float(union_size)

    def merged_tanimoto(self, other):
        """
        Merge Pi-Pi, T-stack, and Hydro
        """
        return self.hydro_overlap(other) / float(self.hydro_footprint() + other.hydro_footprint() - self.hydro_overlap(other))

    def tanimoto_min_norm(self, other):
        """
        Instead of normalizing by the union of bits, normalize by the minimum number of bits
        """
        overlap = self.overlap(other)
        flatten = lambda x: [item for sublist in x for item in sublist]
        total_bits = min(sum(flatten((self.feats.values()))), sum(flatten(other.feats.values())))
        union_size = total_bits
        return overlap / float(union_size)


    def tanimoto_coef_any(self, other):
        """
        Be ambivelent as to the particular type of interaction.
        """
        overlap = 0
        for resi in self.feats:
            if resi in other.feats:
                overlap += 1

        total_bits = min(len(self.feats), len(other.feats))
        union_size = total_bits
        return overlap / float(union_size)

    def merge(self, other):
        """
        Just get the number of overlaps
        """
        overlap = {}
        for resi in self.feats:
            if resi in other.feats:
                overlap[resi] = [i+j for i, j in zip(self.feats[resi], other.feats[resi])]
        return overlap

    def pretty(self):
        out = '\t'.join(['RESI'] + list(self.INTERACTIONS)) + '\n'
        for resi in sorted(map(int, self.feats.keys())):
            out += '\t'.join(map(str, ([resi]+self.feats[resi]))) + '\n'
        return out.strip()

    def num_str(self):
        return':'.join(','.join(map(str, [self.to_common(resi)] + self.feats[resi])) for resi in self.feats)

    def __str__(self):
        return self.pdb + ';' + ':'.join(','.join(map(str, [resi] + self.feats[resi])) for resi in self.feats)

if __name__ == '__main__':
    import sys
    finger_print = FingerPrint.binana_parser(sys.argv[1], sys.argv[2])
    if sys.argv[-1] == 'pretty':
        print finger_print.pretty()
    else:
        print finger_print

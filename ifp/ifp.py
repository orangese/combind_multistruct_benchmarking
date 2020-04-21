import os
import sys

from schrodinger.structure import StructureReader
from schrodinger.structutils.measure import get_shortest_distance

from ifp.hbond import HBond_Container
from ifp.saltbridge import SB_Container
from ifp.hydrophobic import Hydrophobic_Container

nonpolar = {'H': 1.2, 'C':1.7, 'F':1.47, 'Cl':1.75, 'I':1.98, 'Br':1.85}

class AtomGroup:
    def __init__(self, st, st_id):
        self.st = st
        self.st_id = st_id
        self.hacc = [a for a in st.atom if self.valid_acceptor(a)]
        self.hdon = [a for a in st.atom if self.valid_donor(a)]
        self.chrg = [a for a in st.atom if a.formal_charge != 0]
        self.nonpolar = set([a for a in st.atom if self.is_nonpolar(a)])

    def valid_donor(self, atom):
        return (atom.element in ('N', 'O')
                and ('H' in [n.element for n in atom.bonded_atoms])
                and atom.formal_charge >= 0)

    def valid_acceptor(self, atom):
        if atom.element not in ['N', 'O']:
            return False
        return atom.formal_charge <= 0

    def is_nonpolar(self, atom):
        if atom.element == 'H': return False
        return atom.element in nonpolar

class IFP:
    def __init__(self, settings, input_file, output_file, poses):
        self.settings = settings
        self.input_file = input_file
        self.output_file = output_file
        self.poses = poses

        self.protein = {}
        self.ifp = self.fingerprint_pose_viewer()
        self.write_ifp()

    def fingerprint(self, lig_st, prot_st):
        lig = AtomGroup(lig_st, 'lig')

        interactions = {
            'hbond': HBond_Container(lig, [2,3], self.settings),
            'saltbridge': SB_Container(lig, [1], self.settings),
            'hydrophobic': Hydrophobic_Container(lig, [11], self.settings)
        }

        fp = {}
        active_res = []
        for res in prot_st.residue:
            if get_shortest_distance(res.extractStructure(), st2=lig_st, cutoff=8) is not None:
                res_key = '{}:{}({})'.format(res.chain.strip(), res.resnum, res.pdbres.strip())
                if res_key not in self.protein:
                    self.protein[res_key] = AtomGroup(res.extractStructure(), res_key)
                active_res.append(res_key)

        for i_type in interactions:
            for res_key in active_res:
                interactions[i_type].add_residue(res_key, self.protein[res_key])
            
            interactions[i_type].filter_int()
            i_scores = interactions[i_type].score()
            for sc_key, sc in i_scores.items(): 
                fp[sc_key] = sc
        return fp

    def fingerprint_pose_viewer(self):
        with StructureReader(self.input_file) as st:
            prot_st, *poses = list(st)
            poses = poses[:self.poses]

        ifp = []
        for pose in poses:
            ifp += [self.fingerprint(pose, prot_st)]
        return ifp

    def write_ifp(self):
        with open(self.output_file, 'w') as f:
            for i, ifp in enumerate(self.ifp):
                f.write('Pose {}\n'.format(i))
                for sc_key in sorted(ifp.keys()):
                    i,r,ss = sc_key
                    sc = ifp[sc_key]
                    if sc > 0:
                        f.write('{}-{}-{}={}\n'.format(i,r,ss, sc))

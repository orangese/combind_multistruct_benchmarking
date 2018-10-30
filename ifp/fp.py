import os
import sys

from schrodinger.structure import StructureReader
from schrodinger.structutils.measure import get_shortest_distance

from hbond import HBond_Container
from saltbridge import SB_Container
from pipi import PiPi_Container
from hydrophobic import Hydrophobic_Container

nonpolar = {'H': 1.2, 'C':1.7, 'F':1.47, 'Cl':1.75, 'I':1.98, 'Br':1.85}
def valid_donor(atom):
    return (atom.element in ('N', 'O')) and ('H' in [n.element for n in atom.bonded_atoms]) and atom.formal_charge >= 0
def valid_acceptor(atom):
    if atom.element not in ['N', 'O']:
        return False
    return atom.formal_charge <= 0
def is_nonpolar(atom):
    if atom.element == 'H' and [n.element for n in atom.bonded_atoms][0] in nonpolar:
        return True
    return atom.element != 'H' and atom.element in nonpolar
def is_nonpolar2(atom):
    if atom.element == 'H': return False
    return atom.element in nonpolar

class AtomGroup:
    def __init__(self, st, st_id):
        self.st = st
        self.st_id = st_id
        self.hacc = [a for a in st.atom if valid_acceptor(a)]
        self.hdon = [a for a in st.atom if valid_donor(a)]
        self.chrg = [a for a in st.atom if a.formal_charge != 0]
        
        self.cat = [a for a in st.atom if a.formal_charge < 0]
        self.nonpolar1 = set([a for a in st.atom if is_nonpolar(a)])
        self.nonpolar2 = set([a for a in st.atom if is_nonpolar2(a)])

        self.aro = self.fuse_rings(st)

    def fuse_rings(self, st):
        
        aro = [ri for ri in st.ring if ri.isAromatic() or ri.isHeteroaromatic()]

        # Identify pairs of rings sharing atoms
        pairs_to_merge = []
        for i, ring1 in enumerate(aro):
            for j, ring2 in enumerate(aro[i+1:]):
                if any(atom1.xyz == atom2.xyz
                       for atom1 in ring1.atom
                       for atom2 in ring2.atom):
                    pairs_to_merge += [(i, j+i+1)]

        # Merge into groups of disjoint atoms
        to_merge = [set([a]) for a in range(len(aro))]
        for (i, j) in pairs_to_merge:
            group1 = [a for a, group in enumerate(to_merge) if i in group]
            group2 = [a for a, group in enumerate(to_merge) if j in group]

            if group1 == group2: continue
            assert len(group1) == 1
            assert len(group2) == 1
            group1 = group1[0]
            group2 = group2[0]
            if group1 > group2:
                group1, group2 = group2, group1
            to_merge += [to_merge[group1].union(to_merge[group2])]
            to_merge.pop(group2)
            to_merge.pop(group1)

        # Merge structures of joined rings
        fused = []
        for group in to_merge:
            fused += [aro[group.pop()].extractStructure(copy_props=True)]
            for i in group:
                fused[-1].extend(aro[i].extractStructure(copy_props=True))
        return fused


class FP:
    def __init__(self, args):
        self.params = {
            'mode': '',
            'input_file': '',
            'output_file': '',
            'poses': 105,
            'raw': False
        }

        self.set_user_params(args)
        self.protein = {}

        if self.params['mode'] == 'pv':
            self.fp = self.fingerprint_pose_viewer()
        elif self.params['mode'] == 'st':
            self.fp = self.fingerprint_pair()
        self.write_fp()

    def fingerprint(self, lig_st, prot_st, pnum=None):
        lig = AtomGroup(lig_st, 'lig')

        interactions = {
            'hbond': HBond_Container(lig, [2,3]),
            'saltbridge': SB_Container(lig, [1]),
            'pipi': PiPi_Container(lig, [6]),
            'hydrophobic': Hydrophobic_Container(lig, [10,11])
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
            if self.params['raw']:
                i_scores = interactions[i_type].raw()
            else:
                i_scores = interactions[i_type].score()
            for sc_key, sc in i_scores.items(): 
                fp[sc_key] = sc
        return fp

    def fingerprint_pose_viewer(self):
        fp = []
        prot_st = None
        for i, st in enumerate(StructureReader(self.params['input_file'])):
            if i > self.params['poses']: break
            if i == 0:
                prot_st = st
                continue
            fp.append(self.fingerprint(st, prot_st, i))
            
        return fp

    def fingerprint_pair(self):
        pdb = self.params['output_file'].split('_')[0]

        try:
            prot_st = next(StructureReader('../../structures/proteins/{}_prot.mae'.format(pdb)))
        except:
            prot_st = next(StructureReader('../../structures/proteins/{}_prot.mae'.format(os.listdir('../../docking/grids')[0])))

        lig_st = next(StructureReader('../../structures/ligands/{}_lig.mae'.format(pdb)))

        return [self.fingerprint(lig_st, prot_st)]

    def set_user_params(self, args):
        for index in range(len(args)):
            item = args[index]
            if item[0] == '-':
                key = item.replace('-','')
                assert key in self.params, "Key name {} is invalid".format(key)
                if key == 'poses':
                    value = int(args[index+1])
                elif key == 'raw':
                    value = args[index+1] == 'True'
                else:
                    value = args[index+1]
                self.params[key] = value

    def write_fp(self):
        with open(self.params['output_file'], 'w') as f:
            for i, ifp in enumerate(self.fp):
                f.write('Pose {}\n'.format(i))
                for sc_key in sorted(ifp.keys()):#.items():
                    i,r,ss = sc_key
                    sc = ifp[sc_key]
                    if sc >= 0: f.write('{}-{}-{}={}\n'.format(i,r,ss, sc))

if __name__ == '__main__':
    FP(sys.argv[:])


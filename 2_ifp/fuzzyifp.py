import os
import sys
from schrodinger.structure import StructureReader
from schrodinger.structutils.measure import get_shortest_distance

from schrodinger.structutils.rmsd import renumber_conformer

from hbond import HBond_Container
from saltbridge import SB_Container
from pipi import PiPi_Container
from picat import PiCat_Container
from hydrophobic import Hydrophobic_Container
#from vdw import VDW_Container

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
        self.aro = [ri for ri in st.ring if ri.isAromatic() or ri.isHeteroaromatic()]
        self.cat = [a for a in st.atom if a.formal_charge < 0]
        self.nonpolar1 = set([a for a in st.atom if is_nonpolar(a)])
        self.nonpolar2 = set([a for a in st.atom if is_nonpolar2(a)])

class FuzzyIFP:
    def __init__(self, args):
        self.params = {
            'mode': '',
            'input_file': '',
            'output_file': '',
            'verbose':''}
        self.set_user_params(args)
        #print self.params
        #target = open(self.params['verbose'],'a')
        #target.write('\nFingerprinting: ' + self.params['input_file'] + '\n')
        #target.close()

        self.protein = {}

        if self.params['mode'] == 'pv':
            self.fp = self.fingerprint_pose_viewer()
        elif self.params['mode'] == 'st':
            self.fp = self.fingerprint_pair()
        self.write_fp()

    def fingerprint(self, lig_st, prot_st, pnum=None):
        lig = AtomGroup(lig_st, 'lig')
        interactions = {
            #'hal' : HalBond_Container(lig_st, sub_st_map, [1]),
            'hbond': HBond_Container(lig, [2,3]),
            'saltbridge': SB_Container(lig, [0,1,4]),
            'pipi': PiPi_Container(lig, [5,6]),
            'picat': PiCat_Container(lig, [7,8]),
            #'metal': Metal_Container(lig_st, sub_st_map, [9]),
            'hydrophobic': Hydrophobic_Container(lig, [10,11])
            #'vdw': VDW_Container(lig_st)
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
        
        #self.verbose_output(interactions, pnum)

        return fp

    def fingerprint_pose_viewer(self):
        fp = []
        prot_st = None
        
        for i, st in enumerate(StructureReader(self.params['input_file'])):
            if i > 105: break

            if i == 0:
                prot_st = st
                continue
            
            fp.append(self.fingerprint(st, prot_st, i))
            
        return fp

    def fingerprint_pair(self):
        pdb = self.params['output_file'].split('_')[0]

        prot_st = StructureReader('../../structures/proteins/{}_prot.mae'.format(pdb)).next()
        lig_st = StructureReader('../../structures/ligands/{}_lig.mae'.format(pdb)).next()
        
        return [self.fingerprint(lig_st, prot_st)]

    def verbose_output(self, interactions, p_num=None):
        #target = open(self.params['verbose'],'a')

        #if p_num is not None: target.write('\nPose Number '+str(p_num) + '\n')

        #for i_type in interactions:
        #    target.write(str(interactions[i_type])+'\n')

        #for res in sorted(active_site.residues.keys(), key=lambda r:int(r)):
        #    if len(active_site.int_per_res[res]) > 0:
        #        target.write('\nResidue ' + str(res) + '\n')
        #        for i in active_site.int_per_res[res]:
        #            target.write(str(i) + '\n')
        #target.close()
        pass

    def set_user_params(self, args):
        for index in range(len(args)):
            item = args[index]
            if item[0] == '-': # so it's a parameter key value
                key = item.replace('-','')
                assert key in self.params, "Key name {} is invalid".format(key)
                try: value = float(args[index+1])
                except ValueError: value = args[index+1]
                self.params[key] = value

    def write_fp(self):
        with open(self.params['output_file'], 'w') as f:
            for i, ifp in enumerate(self.fp):
                f.write('Pose {}\n'.format(i))
                for sc_key in sorted(ifp.keys()):#.items():
                    i,r,ss = sc_key
                    sc = ifp[sc_key]
                    if sc >= 0.05: f.write('{}-{}-{}={}\n'.format(i,r,ss, sc))
 
    #def __str__(self):
    #    return '\n'.join(';'.join(','.join(map(str, [key,val])) for key, val in fp.items()) for fp in self.fp)

#if __name__ == '__main__':
#    import sys
#    print FuzzyIFP(sys.argv[:])
FuzzyIFP(sys.argv[:])

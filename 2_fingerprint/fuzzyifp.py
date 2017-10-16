#!/share/PI/rondror/software/schrodinger2017-1/run

import os
from schrodinger.structure import StructureReader
from schrodinger.structutils.measure import get_shortest_distance
from schrodinger.structutils.analyze import AslLigandSearcher

from halbond import HalBond_Container
from hbond import HBond_Container
from saltbridge import SB_Container
from pipi import PiPi_Container
from picat import PiCat_Container
from metal import Metal_Container
from hydrophobic import Hydrophobic_Container
from vdw import VDW_Container

class FuzzyIFP:
    def __init__(self, args):
        self.params = {
            'input': '',
            'receptor': '',
            'ligand': '',
            'output_file': '',
            'verbose':''}
        self.set_user_params(args)

        target = open(self.params['verbose'],'a')
        target.write('\nFingerprinting: ' + self.params['receptor'] + '\n')
        target.close()

        if self.params['input'] == 'pose_viewer':
            self.fp = self.fingerprint_pose_viewer()
        else:
            self.fp = [self.fingerprint_pair()]

    def fingerprint(self, lig_st, prot_st, pnum=None):
        interactions = {
            'hal' : HalBond_Container(lig_st),
            'hbond': HBond_Container(lig_st),
            'saltbridge': SB_Container(lig_st),
            'pipi': PiPi_Container(lig_st),
            'picat': PiCat_Container(lig_st),
            'metal': Metal_Container(lig_st),
            'hydrophobic': Hydrophobic_Container(lig_st),
            'vdw': VDW_Container(lig_st)
        }

        num_scores = {'hal':2, 'hbond':2, 'saltbridge':1, 'pipi':1, 'picat':2, 'metal':1,'hydrophobic':2,'vdw':2}

        active_site_res = {}
        for res in prot_st.residue:
            if get_shortest_distance(res.extractStructure(), st2=lig_st)[0] <= 8:
                res_key = '{}:{}({})'.format(res.chain.strip(), res.resnum, res.pdbres.strip())
                active_site_res[res_key] = res.extractStructure()

        fp = {r:[] for r in active_site_res}
        for i_type in ['hal','hbond','saltbridge','pipi','picat','metal','hydrophobic','vdw']: # interactions:
            for r, r_st in active_site_res.items():
                interactions[i_type].add_residue(r, r_st)
            
            interactions[i_type].filter_int()
            i_scores = interactions[i_type].score()
                    
            for r in active_site_res:
                fp[r].extend(i_scores.get(r, [0]*num_scores[i_type] ))
                    
        self.verbose_output(interactions, pnum)

        return fp

    def fingerprint_pose_viewer(self):
        fp = []
        prot_st = None
        for i, st in enumerate(StructureReader(self.params['receptor'])):
            if i > 50: break

            if i == 0:
                prot_st = st
                continue

            fp.append(self.fingerprint(st, prot_st, i))
            
        return fp

    def fingerprint_pair(self):
        prot_st = StructureReader(self.params['receptor']).next()
        
        lig_search = AslLigandSearcher()
        lig_atoms = []
        for l in lig_search.search(prot_st):
            lig_atoms.extend(l.atom_indexes)
        prot_st.deleteAtoms(lig_atoms)

        lig_st = StructureReader(self.params['ligand']).next()
        
        return self.fingerprint(lig_st, prot_st)

    def verbose_output(self, interactions, p_num=None):
        target = open(self.params['verbose'],'a')

        if p_num is not None: target.write('\nPose Number '+str(p_num) + '\n')

        for i_type in interactions:
            target.write(str(interactions[i_type])+'\n')

        #for res in sorted(active_site.residues.keys(), key=lambda r:int(r)):
        #    if len(active_site.int_per_res[res]) > 0:
        #        target.write('\nResidue ' + str(res) + '\n')
        #        for i in active_site.int_per_res[res]:
        #            target.write(str(i) + '\n')
        target.close()

    def set_user_params(self, args):
        for index in range(len(args)):
            item = args[index]
            if item[0] == '-': # so it's a parameter key value
                key = item.replace('-','')
                assert key in self.params, "Key name {} is invalid".format(key)
                try: value = float(args[index+1])
                except ValueError: value = args[index+1]
                self.params[key] = value
    
    def __str__(self):
        return '\n'.join(';'.join(','.join(map(str, [key]+val)) for key, val in fp.items()) for fp in self.fp)

if __name__ == '__main__':
    import sys
    print FuzzyIFP(sys.argv[:])

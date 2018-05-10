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

        #if self.params['input'] == 'pose_viewer':
        if self.params['mode'] == 'pv':
            self.fp = self.fingerprint_pose_viewer()
        elif self.params['mode'] == 'st':
            self.fp = self.fingerprint_pair()
        self.write_fp()
        #else:
        #    self.fp = [self.fingerprint_pair()]

    def fingerprint(self, lig_st, prot_st, sub_st_map, pnum=None):
        interactions = {
            #'hal' : HalBond_Container(lig_st, sub_st_map, [1]),
            'hbond': HBond_Container(lig_st, sub_st_map, [2,3]),
            'saltbridge': SB_Container(lig_st, sub_st_map, [4]),
            'pipi': PiPi_Container(lig_st, sub_st_map, [5,6]),
            'picat': PiCat_Container(lig_st, sub_st_map, [7,8]),
            #'metal': Metal_Container(lig_st, sub_st_map, [9]),
            'hydrophobic': Hydrophobic_Container(lig_st, sub_st_map,[10,11])
            #'vdw': VDW_Container(lig_st)
        }

        #num_scores = {'hal':2, 'hbond':2, 'saltbridge':1, 'pipi':2, 'picat':2, 'metal':1,'hydrophobic':2}#,'vdw':2}

        active_site_res = {}
        for res in prot_st.residue:
            if get_shortest_distance(res.extractStructure(), st2=lig_st)[0] <= 8:
                res_key = '{}:{}({})'.format(res.chain.strip(), res.resnum, res.pdbres.strip())
                active_site_res[res_key] = res.extractStructure()

        fp = {} # {r:[] for r in active_site_res}
        for i_type in interactions:# ['hal','hbond','saltbridge','pipi','picat','metal','hydrophobic']:#, 'other']:#,'vdw']: # interactions:
            if i_type != 'other':
                for r, r_st in active_site_res.items():
                    interactions[i_type].add_residue(r, r_st)
            
                interactions[i_type].filter_int()
                i_scores = interactions[i_type].score()
                for sc_key, sc in i_scores.items(): fp[sc_key] = sc     
            #for r in active_site_res:
                #if i_type == 'other':
                #    fp[r].extend([fp[r][0] + fp[r][2], fp[r][1] + fp[r][3]])
                #else:
                
            #    fp[r].extend(i_scores.get(r, [0]*num_scores[i_type] ))
        #self.verbose_output(interactions, pnum)

        return fp

    def fingerprint_pose_viewer(self):
        fp = []
        prot_st = None

        #os.chdir(self.params['data_dir'])

        #ref_ligand = StructureReader(self.params['ligand']).next()
        ref_ligand = None
        sub_st_map = {}

        for i, st in enumerate(StructureReader(self.params['input_file'])):
            if i > 105: break

            if i == 0:
                prot_st = st
                continue

            renumbered_st = st
            if ref_ligand is not None: 
                renumbered_st = renumber_conformer(ref_ligand, st, use_symmetry=True)

            fp.append(self.fingerprint(renumbered_st, prot_st, sub_st_map, i))
            
        return fp

    def fingerprint_pair(self):
        pdb = self.params['output_file'].split('_')[0]

        prot_st = StructureReader('../../structures/proteins/{}_prot.mae'.format(pdb)).next()
        lig_st = StructureReader('../../structures/ligands/{}_lig.mae'.format(pdb)).next()
        
        #lig_search = AslLigandSearcher()
        #lig_atoms = []
        #for l in lig_search.search(prot_st):
        #    lig_atoms.extend(l.atom_indexes)
        #prot_st.deleteAtoms(lig_atoms)

        #lig_st = StructureReader(self.params['ligand']).next()
        
        return [self.fingerprint(lig_st, prot_st, {})]

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

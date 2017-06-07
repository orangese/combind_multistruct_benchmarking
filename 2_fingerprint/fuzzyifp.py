#!/share/PI/rondror/software/schrodinger2017-1/run
# Change this to /share/PI/rondror/software/schrodinger2016-1/run? Testing new stuff out
from interactions import Interactions
import os
from schrodinger import structure # Make sure you are using schrodinger python distribution!
from schrodinger.structutils.minimize import minimize_structure
from schrodinger.structutils.analyze import AslLigandSearcher
from receptor import Receptor
from ligand import Ligand
import numpy as np
from tqdm import tqdm
from multiprocessing import Pool

class FuzzyIFP:
    def __init__(self, args):
        self.params = {
            'electrostatic_dist_cutoff': 4.0,
            'hydrophobic_dist_cutoff': 4.5,
            'hydrogen_bond_dist_cutoff': 4.0,
            'hydrogen_bond_angle_cutoff': 40.0,
            'pi_padding_dist': 0.75,
            'pi_pi_interacting_dist_cutoff': 7.5,
            'pi_stacking_angle_tolerance':30.0,
            'T_stacking_angle_tolerance': 30.0,
            'T_stacking_closest_dist_cutoff': 5.0,
            'cation_pi_dist_cutoff': 6.0,
            'input': '',
            'receptor': '',
            'ligand': '',
            'output_file': ''}
        self.set_user_params(args)
        self.struct = structure.StructureReader(self.params['receptor'])
        self.receptor = Receptor()
        st = self.struct.next()
        minimize_structure(st, max_steps = 0) #Assigns partial charges
        
        self.receptor.load_mae(st)

        if self.params['input'] == 'pose_viewer':
            self.fp = self.fingerprint_pose_viewer()
        elif self.params['input'] == 'complex':
            self.fp = [self.fingerprint_complex()]
        else:
            self.fp = [self.fingerprint_pair()]
        
    def get_fingerprint(self, values):
        index, mergedReceptor, lig = values
        return (index, Interactions(mergedReceptor, lig, self.params).get_fp())
    
    def fingerprint_pose_viewer(self):
        fp = []
        
        for pose_num, st in enumerate(self.struct):
            if pose_num > 50: #Only fingerprint poses 0->50, right now this is just a time saving measure
                continue

            mergedReceptorSt = self.receptor.st.merge(st)
            minimize_structure(mergedReceptorSt, max_steps = 0)
            mergedReceptor = Receptor()
            mergedReceptor.load_mae(mergedReceptorSt)
            
            lig = Ligand()
            lig.load_mae(st)
            fp.append(Interactions(mergedReceptor, lig, self.params).get_fp())
            
        return fp


    def fingerprint_pair(self):
        lig = Ligand()
        st = structure.StructureReader(self.params['ligand']).next()
        minimize_structure(st, max_steps = 0)
        lig.load_mae(st) #Assigns partial charges
        return Interactions(self.receptor, lig, self.params).get_fp()

    def fingerprint_complex(self):
        lig = self.receptor.export_ligand(self.params['ligand'])
        return Interactions(self.receptor, lig, self.params).get_fp()
                     
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
        return '\n'.join(':'.join(','.join(map(str, [key]+val)) for key, val in fp.items()) for fp in self.fp)

if __name__ == '__main__':
    import sys
    print FuzzyIFP(sys.argv[:])

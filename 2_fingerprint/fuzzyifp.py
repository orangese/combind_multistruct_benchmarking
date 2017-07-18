#!/share/PI/rondror/software/schrodinger2017-1/run
#from interactions import Interactions
import os
from schrodinger import structure # Make sure you are using schrodinger python distribution!
from schrodinger.structutils.minimize import minimize_structure
from receptor import Receptor
from point_atom import Atom
from ligand import Ligand

class FuzzyIFP:
    def __init__(self, args):
        self.params = {
            'active_site_cutoff': '5', # angstroms
            'input': '',
            'receptor': '',
            'ligand': '',
            'output_file': '',
            'verbose':'',
            'receptor_name':''}
        self.set_user_params(args)
        
        self.struct = structure.StructureReader(self.params['receptor'])
        self.receptor = Receptor()
        st = self.struct.next()
        minimize_structure(st, max_steps = 0) #Assigns partial charges
        self.receptor.load_mae(st)

        target = open(self.params['verbose'],'a')
        target.write('\nFingerprinting: ' + self.params['receptor'] + '\n')
        target.close()

        if self.params['input'] == 'pose_viewer':
            self.fp = self.fingerprint_pose_viewer()
        else:
            self.fp = [self.fingerprint_pair()]
        
    def get_active_site(self, ligand, receptor):
        active_site = Receptor()
        #Only add residues that are within RESIDUE_THRESH of the ligand center coordinate
        for residue in receptor.residues.values():
            added = False
            for ligand_atom in ligand.atoms:
                for receptor_atom in residue.atoms:
                    if ligand_atom.coordinates.dist_to(receptor_atom.coordinates) < self.params['active_site_cutoff']:
                        assert not added
                        active_site.add_residue(residue.copy_of())
                        added = True
                        break
                if added: break
         
        active_site.assign()
        
        #Calculate Atom ID -> Res, based on the merged structure
        active_site.st = receptor.st
        atomIDToRes = {}
        for atom in active_site.st.atom:
            cur_atom = Atom()
            cur_atom.from_schrod(atom)
            atomIDToRes[cur_atom.atom_id] = cur_atom.residue_id
        active_site.atomIDToRes = atomIDToRes
        active_site.seperate_atoms(receptor, ligand)
        return active_site

    def fingerprint_pose_viewer(self):
        fp = []
 
        for pose_num, st in enumerate(self.struct):
            mergedReceptorSt = self.receptor.st.merge(st)
            minimize_structure(mergedReceptorSt, max_steps = 0)
            mergedReceptor = Receptor()
            mergedReceptor.load_mae(mergedReceptorSt)
            
            lig = Ligand()
            lig.load_mae(st)

            active_site = self.get_active_site(lig, mergedReceptor)

            fp.append(active_site.fingerprint(lig))

            self.verbose_output(active_site, p_num=pose_num)
            
        return fp

    def fingerprint_pair(self):
        lig = Ligand()
        st = structure.StructureReader(self.params['ligand']).next()
        minimize_structure(st, max_steps = 0)
        lig.load_mae(st) #Assigns partial charges

        active_site = self.get_active_site(lig, self.receptor)
        
        fp = active_site.fingerprint(lig)

        self.verbose_output(active_site)
        self.print_size(lig)
        return fp

    def print_size(self, ligand):
        weight = ligand.molecular_weight()
        diameter = ligand.diameter()

        target = open('/scratch/PI/rondror/docking_data/'+self.params['receptor_name']+'/lig_sizes.txt', 'a')
        target.write(self.params['ligand'] + ',' + str(weight) + ',' + str(diameter) + '\n')
        target.close()

    def verbose_output(self, active_site, p_num=None):
        target = open(self.params['verbose'],'a')

        if p_num is not None: target.write('\nPose Number '+str(p_num) + '\n')

        for res in sorted(active_site.residues.keys()):
            if len(active_site.residues[res].all_interactions) > 0:
                target.write('\nResidue ' + str(res) + '\n')
                for i in active_site.residues[res].all_interactions:
                    target.write(str(i) + '\n')
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
        return '\n'.join(':'.join(','.join(map(str, [key]+val)) for key, val in fp.items()) for fp in self.fp)

if __name__ == '__main__':
    import sys
    print FuzzyIFP(sys.argv[:])

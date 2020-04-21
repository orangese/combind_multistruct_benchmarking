import os
from glob import glob

from utils import grouper
from mcss.mcss import MCSS

class MCSSController:
    """
    Controls computation and loading of MCSS features for a given target.

     There are two key phases of the computation:
        (1) Identification of maximum common substructure(s)
        (2) Computation of RMSDs between the substructures in
            docking results.

    The role of this class is to decide which ligand pairs should have MCSSs
    computed, and, for the chosen few, how far along they are in the process.
    This class also sets all file and path names, which are then passed to
    mcss.py for the actual computation.

    This class collates outputs from individual MCSS computations
    into a single file. This is important, as otherwise you end up making
    on order of # PDB * # CHEMBL files!

    lm: LigandManager

    MCSSs: {ligand-pair-string: MCSS instance, ...}
    pdb: PDB x-docked pdb ligands
    docked: set(string), Set of ligands with docking results
    no_mcss: set(string), Ligands with no init MCSS
    no_rmsd: set(string), Ligands with no MCSS RMSD
    no_mcss_small: set(string), Ligands that are small enough that the
        normal lower bound on MCSS size could not be met while a significant
        overlap exists. Mostly for use in making a plot of ligand similarity
        versus docking performance.
    """

    INIT_GROUP_SIZE = 500
    RMSD_GROUP_SIZE = 5

    QUEUE = 'rondror'
    
    TEMPLATE = ('#!/bin/bash\n'
                '#SBATCH -p {}\n'
                '#SBATCH --tasks=1\n'
                '#SBATCH -t 2:00:00\n'
                '{}'
                'wait\n'
                ).format(QUEUE, '{}')


    def __init__(self, lm):
        self.lm = lm
        self.pdb = lm.get_xdocked_ligands(lm.params['n_ligs'])
        self.docked  = set(lm.docked(self.pdb))

        self.MCSSs = {}
        self.no_mcss = set([])
        self.no_rmsd = set([])

        # File paths. Remaining '{}' is ligand pair name.
        self.root = lm.path('MCSS')
        self.atom_types = '{}/mcss/custom_types/{}.typ'.format(lm.path('CODE'),
                                                               lm.params['mcss_version'])
        self.mcss_file = "{}/mcss.csv".format(self.root)
        self.init_file = "{}/{}.init.csv".format(self.root, '{}')

        self.rmsd_file = "{}/{}-{}-{}.csv".format(self.root, '{}',
                                                  self.lm.st, lm.params['docking_version'])

        self.init_command, self.rmsd_command = self._construct_commands(
                                                    lm.params['max_poses'])

    def _construct_commands(self, max_poses):
        """
        Construct commands to run MCSS computation.
        
        Commands contain remaining substitutions to handle ligand specific info:
            '{0:},{1:}' ligand1,ligand2, '{2:},{3:}' poses1,poses2 '{4:}' str(MCSS).
        """
        
        init_command =  'run {0:}/main.py mcss INIT '.format(self.lm.path('CODE'))
        init_command += '{0:} {1:} {2:} {3:} ' # ligand names
        init_command += self.init_file.format('{0:}-{1:}') + ' '
        init_command += self.atom_types + ' '
        init_command += '\n'
        
        rmsd_command = 'run {0:}/main.py mcss RMSD '.format(self.lm.path('CODE'))
        rmsd_command += '{0:} {1:} {2:} {3:} '
        rmsd_command += self.init_file.format('{0:}-{1:}') + ' '
        rmsd_command += self.atom_types + ' '
        
        rmsd_command += self.rmsd_file.format('{0:}-{1:}')
        rmsd_command += " {} ".format(max_poses)
        rmsd_command += ' "{4:}"\n'

        return init_command, rmsd_command

    # Primary public methods
    def get_rmsd(self, l1, l2, pose1, pose2):
        """
        Gets the rmsd for a pair of poses.

        l1, l2: string, ligand names
        pose1, pose2: int, pose numbers for which to get rmsd.
        """
        if l1 > l2:
            l1, l2 = l2, l1
            pose1, pose2 = pose2, pose1
        mcss = self.MCSSs["{}-{}".format(l1, l2)]
        if not self.is_valid(mcss): return None
        if (pose1, pose2) not in mcss.rmsds: return None
        return mcss.rmsds[(pose1, pose2)]

    def load_rmsds(self, ligands, max_poses):
        """
        Load rmsds into memory. Will error if any requested rmsds have
        not been computed.

        ligands: iterator, set of ligands for which to load rmsds
        max_poses: int, maximum number of poses for which to load rmsds
        """
        self.load_mcss()
        ligands = list(set(ligands))
        for i, _l1 in enumerate(ligands):
            for _l2 in ligands[i+1:]:
                if _l1 > _l2:
                    l1, l2 = _l2, _l1
                else:
                    l1, l2 = _l1, _l2
                name = "{}-{}".format(l1, l2)
                assert name in self.MCSSs, name
                if self.is_valid(self.MCSSs[name]):
                    poseviewer_paths = {l1: self.lm.path('DOCK_PV', {'ligand': l1}),
                                        l2: self.lm.path('DOCK_PV', {'ligand': l2})}
                    self.MCSSs[name].load_rmsds(self.rmsd_file.format(name),
                                                poseviewer_paths, max_poses)

    def get_mcss_and_ligand_sizes(self, l1, l2):
        """
        Returns the fraction of the smaller of the ligands
        that composes the MCSS for l1 and l2.

        l1, l2: string, ligand names
        """
        if l1 > l2: l1, l2 = l2, l1
        mcss = self.MCSSs["{}-{}".format(l1, l2)]
        return mcss.n_mcss_atoms, mcss.n_l1_atoms, mcss.n_l2_atoms

    def get_mcss_size(self, l1, l2, compute=False):
        """
        Returns the fraction of the smaller of the ligands
        that composes the MCSS for l1 and l2.

        l1, l2: string, ligand names
        """
        if l1 > l2: l1, l2 = l2, l1
        try:
            mcss = self.MCSSs["{}-{}".format(l1, l2)]
        except KeyError:
            self._execute_init_on_the_fly(l1, l2)
            mcss = self.MCSSs["{}-{}".format(l1, l2)]
        f = self.lm.params['mcss_func']
        return mcss.n_mcss_atoms / float(f(mcss.n_l1_atoms, mcss.n_l2_atoms))

    def is_valid(self, mcss):
        """
        Decides if the mcss is large enough to include in calculations.
        
        (It would be problematic if small MCSSs were taken into account
        because the score is based solely on RMSD).
        """
        f = self.lm.params['mcss_func']
        r = self.lm.params['mcss_rel_min']
        rel_thresh = r * f(mcss.n_l1_atoms, mcss.n_l2_atoms)
        
        abs_thresh = self.lm.params['mcss_abs_min']

        return mcss.n_mcss_atoms > max(rel_thresh, abs_thresh)

    def load_mcss(self, temp_init_files=None):
        """
        Read MCSSs both from directories and from collated file. Does not
        load rmsds! Do this with load_rmsd method.

        temp_init_files: string, file names of temporary mcss files.
        """
        self.MCSSs = {}
        if os.path.exists(self.mcss_file):
            with open(self.mcss_file) as fp:
                for line in fp:
                    try:
                        mcss = MCSS.from_string(line)
                    except ValueError:
                        assert False, line
                    assert mcss.name not in self.MCSSs, mcss.name
                    self.MCSSs[mcss.name] = mcss

        if temp_init_files is None:
            temp_init_files = glob('{}/*.init.csv'.format(self.root))
        
        for temp in temp_init_files:
            self._load_temp(temp)

    def _load_temp(self, temp):
        with open(temp) as fp:
            try:
                mcss = MCSS.from_string(fp.readline().strip())
            except ValueError:
                assert False, temp
            if mcss is not None:
                self.MCSSs[mcss.name] = mcss

    # All of below are relevant for computation only
    def compute_mcss(self):
        """
        Compute unfinished MCSS features. See above class description for more detail.
        """
        previous_cwd = os.getcwd()
        os.system('mkdir -p {}'.format(self.root))
        os.chdir(self.root)

        print("{} ligands".format(len(self.pdb)))

        self._collate_mcss()
        self._add_pdb_to_pdb()
        self._execute()
        os.chdir(previous_cwd)
    
    def _collate_mcss(self):
        """
        Collate MCSS init results into a single file. This also loads
        the data and should be used as the load function during the
        processing step.
        """
        # Get temp_init_files before loading, so that we only
        # delete temp_init_files that have been loaded.
        self.MCSSs = {}
        temp_init_files = glob('{}/*.init.csv'.format(self.root))
        self.load_mcss(temp_init_files)

        if not temp_init_files: return

        print("Collating {} init files".format(len(temp_init_files)))

        # Write to temp and then overwrite original file
        # so that we don't lose anything if job crashes mid-run.
        with open('mcss_temp.csv', 'w') as fp:
            for mcss in self.MCSSs.values():
                fp.write(str(mcss)+'\n')
        
        os.system('mv mcss_temp.csv {}'.format(self.mcss_file))
        
        for temp in temp_init_files:
            os.system('rm {}'.format(temp))

    # Methods to add ligand pairs
    def _add_pdb_to_pdb(self):
        """
        Add all pdb to pdb ligand pairs
        """
        compute_rmsds = True
        for i, l1 in enumerate(self.pdb):
            for l2 in self.pdb[i+1:]:
                self._add(l1, l2, compute_rmsds)

    def _add(self, l1, l2, compute_rmsd):
        """
        Check if pair has been computed, if not add it.

        l1: string, ligand name
        l2: string, ligand name
        rmsd: bool, compute RMSD if true
        """
        if l1 > l2: l1, l2 = l2, l1
        name = '{}-{}'.format(l1, l2)
        if name not in self.MCSSs:
            self.no_mcss.add((l1,l2))
        elif (    compute_rmsd
              and (l1 in self.docked)
              and (l2 in self.docked)
              and self.is_valid(self.MCSSs[name])
              and not os.path.exists(self.rmsd_file.format(name))):
            self.no_rmsd.add((l1, l2))

    # Methods to execute computation
    def _execute(self):
        """
        Execute all incomplete computations.
        """
        if self.no_mcss:
            print(len(self.no_mcss), 'mcss init pairs left')
            self._execute_init()
        if self.no_rmsd:
            print(len(self.no_rmsd), 'mcss rmsd pairs left')
            self._execute_rmsd()

    def _execute_init(self, inline=False):
        for i, pairs in enumerate(grouper(self.INIT_GROUP_SIZE, self.no_mcss)):
            script = 'init{}.sh'.format(i)
            contents = 'export SCHRODINGER_CANVAS_MAX_MEM=1e+12\n'
            for l1, l2 in pairs:
                poses1 = self.lm.path('LIGANDS', {'ligand': l1.replace('_crystal', '')})
                poses2 = self.lm.path('LIGANDS', {'ligand': l2.replace('_crystal', '')})
                contents += self.init_command.format(l1, l2, poses1, poses2)
            with open(script, 'w') as f:
                f.write(self.TEMPLATE.format(contents))

            if inline:
                os.system('sh {} 2>&1 > /dev/null'.format(script))
            else:
                os.system('sbatch {}'.format(script))
        self.no_mcss = set([])

    def _execute_rmsd(self):
        for i, pairs in enumerate(grouper(self.RMSD_GROUP_SIZE, self.no_rmsd)):
            script = 'rmsd{}.sh'.format(i)
            contents = ''
            for l1, l2 in pairs:
                mcss = self.MCSSs["{}-{}".format(l1, l2)]
                poses1 = (self.lm.path('LIGANDS', {'ligand': l1.replace('_crystal', '')})
                          if 'crystal' in l1 else
                         self.lm.path('DOCK_PV', {'ligand': l1}))
                poses2 = (self.lm.path('LIGANDS', {'ligand': l2.replace('_crystal', '')})
                          if 'crystal' in l2 else
                          self.lm.path('DOCK_PV', {'ligand': l2}))
                contents += self.rmsd_command.format(l1, l2, poses1, poses2, str(mcss))
            with open(script, 'w') as f:
                f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script))
        self.no_rmsd = set([])

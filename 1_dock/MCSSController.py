import os
import sys
from grouper import grouper

sys.path.append('../2_ifp')
from mcss_main import MCSS

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
    2_ifp/mcss_main.py for the actual computation.

    Additionally, this class collates outputs from individual MCSS computations
    into a single file. This is important, as otherwise you end up making
    on order of # PDB * # CHEMBL files!

    lm: LigandManager
    max_pdb: int, maximum number of pdb ligands to consider
    max_poses: int, maximum number of poses to consider for RMSD calculations

    MCSSs: {ligand-pair-string: MCSS instance, ...}
    docked: set(string), Set of ligands with docking results
    no_mcss: set(string), Ligands with no init MCSS
    no_rmsd: set(string), Ligands with no MCSS RMSD
    """

    INIT_GROUP_SIZE = 100
    RMSD_GROUP_SIZE = 10

    QUEUE = 'owners'
    
    TEMPLATE = ('#!/bin/bash\n'
                '#SBATCH -p {}\n'
                '#SBATCH --tasks=1\n'
                '#SBATCH -t 8:00:00\n'
                'ml load chemistry\n'
                'ml load schrodinger\n'
                '{}'
                'wait\n'
                ).format(QUEUE, '{}')

    def __init__(self, lm, max_pdb=20, max_poses = 100):
        self.lm = lm
        self.max_pdb = max_pdb

        self.MCSSs = {}
        
        self.pdb = self.lm.docked(self.lm.pdb)[:self.max_pdb]
        self.docked  = set(lm.docked(lm.pdb+lm.chembl()))
        
        # Incomplete pairs
        self.no_mcss = set([])
        self.no_rmsd = set([])

        # File paths. Remaining '{}' is ligand pair name.
        self.root = "{}/{}/mcss/{}".format(self.lm.sp['data'], self.lm.prot, self.lm.sp['mcss'])
        self.mcss_file = "{}/{}_mcss.csv".format(self.root, self.lm.prot)
        self.rmsd_file = "{}/{}-{}-{}.csv".format(self.root, '{}',
                                                     self.lm.st, self.lm.sp['docking'])
        
        self.init_command, self.rmsd_command = self._construct_commands(max_poses)

    def _construct_commands(self, max_poses):
        """
        Construct commands to run MCSS computation.
        
        Commands contain remaining substitutions to handle ligand specific info:
            '{0:}' ligand1, '{1:}' ligand2, '{2:}' str(MCSS).
        """
        lig_template = '{0:}/{1:}/ligands/prepared_ligands/{2:}/{2:}.mae '
        pv_template = '{0:}/{1:}/docking/{2:}/{3:}-to-{4:}/{3:}-to-{4:}_pv.maegz '
        
        init_command  = 'mkdir -p {0:}-{1:}; cd {0:}-{1:}; '
        init_command += '$SCHRODINGER/run {0:}/2_ifp/mcss_main.py INIT '.format(self.lm.sp['code'])
        init_command += '{0:} {1:} ' # ligand names
        init_command += lig_template.format(self.lm.sp['data'], self.lm.prot, '{0:}') # l1_path
        init_command += lig_template.format(self.lm.sp['data'], self.lm.prot, '{1:}') # l2_path
        init_command += '{}/2_ifp/custom_types/{}.typ'.format(self.lm.sp['code'],
                                                              self.lm.sp['mcss']) # Atom typing
        init_command += ' ; cd ..\n'
        
        rmsd_command = '$SCHRODINGER/run {0:}/2_ifp/mcss_main.py RMSD '.format(self.lm.sp['code'])
        rmsd_command += '{0:} {1:} ' # ligand names
        rmsd_command += pv_template.format(self.lm.sp['data'], self.lm.prot,
                                           self.lm.sp['docking'], '{0:}', self.lm.st) # pv1_path
        rmsd_command += pv_template.format(self.lm.sp['data'], self.lm.prot,
                                           self.lm.sp['docking'], '{1:}', self.lm.st) # pv2_path
        rmsd_command += self.rmsd_file.format('{0:}-{1:}')
        rmsd_command += " {} ".format(max_poses)
        rmsd_command += ' "{2:}"\n'
        return init_command, rmsd_command

    # Primary public methods
    def get_rmsd(self, l1, l2, pose1, pose2):
        """
        Gets the rmsd for a pair of poses.

        l1, l2: string, ligand names
        pose1, pose2: int, pose numbers for which to get rmsd.
        """
        if l1 > l2: l1, l2 = l2, l1
        return self.MCSSs["{}-{}".format(l1, l2)].rmsds[(p1, p2)]

    def load_rmsds(self, ligands, max_poses):
        """
        Load rmsds into memory. Will error if any requested rmsds have
        not been computed.

        ligands: iterator, set of ligands for which to load rmsds
        max_poses: int, maximum number of poses for which to load rmsds
        """
        self.load_mcss()
        ligands = list(set(ligands))
        for i, ligand1 in enuemrate(ligands):
            for ligand2 in ligands[i+1:]:
                if l1 > l2: l1, l2 = l2, l1
                name = "{}-{}".format(l1, l2)
                assert name in self.MCSSs
                if self.MCSSs[name].is_valid():
                    self.MCSSs[name].load_rmsds(self.rmsd_file.format(name), max_poses)

    def sort_by_mcss(self, query, ligands):
        """
        Sort ligands by size of maximum common substructure with query.
        Will give an error if any pairs have not been computed.

        * Must first call load_mcss *

        query: string, ligand name
        ligands: list(ligands), set of ligands to be sorted
        """
        def size(ligand):
            if ligand < query:
                name = '{}-{}'.format(ligand, query)
            else:
                name = '{}-{}'.format(query, ligand)
            return self.MCSSs[name].n_mcss_atoms
        return sorted(ligands, key=size, reverse=True)

    def load_mcss(self, temp_directories = None):
        """
        Read MCSSs both from directories and from collated file. Does not
        load rmsds! Do this with load_rmsd method.

        temp_directories: string, directory names of temporary mcss files.
        """
        self.MCSSs = {}
        if os.path.exists(self.mcss_file):
            with open(self.mcss_file) as fp:
                for line in fp:
                   mcss = MCSS.from_string(line)
                   assert mcss.name not in self.MCSSs, mcss.name
                   self.MCSSs[mcss.name] = mcss

        if temp_directories is None:
            temp_directories = self._get_temp_directories()
        for temp in temp_directories:
            mcss = MCSS.from_canvasMCSS(temp)
            if mcss is not None:
                self.MCSSs[mcss.name] = mcss

    # All of below are relevant for computation only
    def collate_mcss(self):
        """
        Collate MCSS init results into a single file. This also loads
        the data and should be used as the load function during the
        processing step.
        """
        # Get temp_directories before loading, so that we only
        # delete temp directories that have been loaded.
        self.MCSSs = {}
        temp_directories = self._get_temp_directories()
        self.load_mcss(temp_directories)

        if not temp_directories: return

        # Write to temp and then overwrite original file
        # so that we don't lose anything if job crashes mid-run.
        # (okay, okay, I guess we could lose things during the move...)
        with open('mcss_temp.csv', 'w') as fp:
            for mcss in self.MCSSs.values():
                fp.write(str(mcss)+'\n')
        
        os.system('mv mcss_temp.csv {}'.format(self.mcss_file))
        for temp in temp_directories:
            os.system('rm -rf {}'.format(temp))

    def _get_temp_directories(self):
        """
        Return list of directories with completed MCSS init.
        """
        return [fname for fname in os.listdir(self.root)
                if os.path.isdir("{}/{}".format(self.root, fname))]

    # Methods to add ligand pairs
    def add_pdb_to_pdb(self, compute_rmsds):
        """
        Add all pdb to pdb ligand pairs
        """
        for i, l1 in enumerate(self.pdb):
            for l2 in self.pdb[i+1:]:
                self._add(l1, l2, compute_rmsds)

    def add_pdb_to_allchembl(self):
        """
        Add all pdb - chembl ligand pairs.
        Most of these are only used for deciding which chembl ligands to use
        so don't compute rmsd for all of them.
        """
        compute_rmsds = False
        for l1 in self.pdb:
            for l2 in self.lm.chembl():
                self._add(l1,l2, compute_rmsds)

    def add_pick_helpers(self, pick_helpers, compute_rmsds):
        """
        Adds pairs of pdb to chembl ligands and chembl to chembl ligands
        that will be jointly used in a score computation as specified
        by pick_helpers.

        Chembl - Chembl pairs that are specified by "chembl"
        pick_helpers: {pick_helpers_filename: {query_pdb: [chembl, ...]}}
        """
        for fname, queries in pick_helpers.items():
            for query, chembl_ligands in queries.items():
                for i, l1 in enumerate(chembl_ligands):
                    assert l1 != '' # switch to continue if this happens
                    self._add(query, l1, compute_rmsds)
                    for l2 in chembl_ligands[i+1:]:
                        assert l2 != '' # switch to continue if this happens
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
              and l1 in self.docked and l2 in self.docked
              and self.MCSSs[name].is_valid()
              and not os.path.exists(self.rmsd_file.format(name))):
            self.no_rmsd.add((l1, l2))

    # Methods to execute computation
    def execute(self):
        """
        Execute all incomplete computations.
        """
        if len(self.no_mcss) > 0:
            print(len(self.no_mcss), 'mcss init pairs left')
            self._execute_init()
        if len(self.no_rmsd) > 0:
            print(len(self.no_rmsd), 'mcss rmsd pairs left')
            self._execute_rmsd()

    def _execute_init(self):
        for i, pairs in enumerate(grouper(self.INIT_GROUP_SIZE, self.no_mcss)):
            script = 'init{}.sh'.format(i)
            contents = 'export SCHRODINGER_CANVAS_MAX_MEM=1e+12\n'
            for l1, l2 in pairs:
                contents += self.init_command.format(l1, l2)
            with open(script, 'w') as f: f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script))

    def _execute_rmsd(self):
        for i, pairs in enumerate(grouper(self.RMSD_GROUP_SIZE, self.no_rmsd)):
            script = 'rmsd{}.sh'.format(i)
            contents = ''
            for l1, l2 in pairs:
                contents += self.rmsd_command.format(l1, l2, str(self.MCSSs["{}-{}".format(l1, l2)]))
            with open(script, 'w') as f: f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script))


def compute_mcss(lm, pick_helpers={}, max_pdb=20, max_poses = 100, compute_rmsds = True):
    """
    Compute unfinished MCSS features. See above class description for more detail.

    * This should be called from the root of a protein directory *

    lm: LigandManager instance
    pick_helpers: {pick_helpers_filename: {query_pdb: [chembl, ...]}}
    max_pdb: int, maximum number of pdb ligands to consider
    max_poses: int, maximum number of poses for which to compute rmsds
    """
    os.system('mkdir -p mcss/{}'.format(lm.sp['mcss']))
    os.chdir('mcss/{}'.format(lm.sp['mcss']))
    
    controller = MCSSController(lm, max_pdb, max_poses)
    controller.collate_mcss()
    controller.add_pdb_to_pdb(compute_rmsds)
    controller.add_pdb_to_allchembl()
    controller.add_pick_helpers(pick_helpers, compute_rmsds)
    controller.execute()
    os.chdir('../..')

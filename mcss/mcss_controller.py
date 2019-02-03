import os
from glob import glob

from grouper import grouper
from shared_paths import shared_paths
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
    chembl: all chembl ligands
    docked: set(string), Set of ligands with docking results
    no_mcss: set(string), Ligands with no init MCSS
    no_rmsd: set(string), Ligands with no MCSS RMSD
    no_mcss_small: set(string), Ligands that are small enough that the
        normal lower bound on MCSS size could not be met while a significant
        overlap exists. Mostly for use in making a plot of ligand similarity
        versus docking performance.
    """

    INIT_GROUP_SIZE = 500
    RMSD_GROUP_SIZE = 15

    QUEUE = 'owners'
    
    TEMPLATE = ('#!/bin/bash\n'
                '#SBATCH -p {}\n'
                '#SBATCH --tasks=1\n'
                '#SBATCH -t 2:00:00\n'
                '{}'
                'wait\n'
                ).format(QUEUE, '{}')


    def __init__(self, lm):
        self.st = lm.st
        self.pdb = lm.get_xdocked_ligands(shared_paths['stats']['n_ligs'])
        self.chembl = lm.chembl()
        self.docked  = set(lm.docked(self.pdb+self.chembl))

        self.MCSSs = {}
        self.no_mcss = set([])
        self.no_rmsd = set([])
        self.no_mcss_small = set([])

        # File paths. Remaining '{}' is ligand pair name.
        self.root = "{}/mcss/{}".format(lm.root, shared_paths['mcss'])
        self.atom_types = '{}/mcss/custom_types/{}.typ'.format(shared_paths['code'],
                                                               shared_paths['mcss'])
        self.mcss_file = "{}/mcss.csv".format(self.root)
        self.init_file = "{}/{}.init.csv".format(self.root, '{}')

        self.rmsd_file = "{}/{}-{}-{}.csv".format(self.root, '{}',
                                                  self.st, shared_paths['docking'])

        self.lig_template = '{0:}/ligands/prepared_ligands/{1:}/{1:}.mae'.format(lm.root, '{0:}')
        self.crystal_template = '{0:}/structures/ligands/{1:}.mae'.format(lm.root, '{0:}')
        self.pv_template = '{0:}/docking/{1:}/{2:}-to-{3:}/{2:}-to-{3:}_pv.maegz'.format(
                                lm.root, shared_paths['docking'], '{0:}', lm.st)
        self.init_command, self.rmsd_command = self._construct_commands(
                                                    shared_paths['stats']['max_poses'])

        assert not any('CHEMBL' in ligand for ligand in self.pdb)

    def _construct_commands(self, max_poses):
        """
        Construct commands to run MCSS computation.
        
        Commands contain remaining substitutions to handle ligand specific info:
            '{0:},{1:}' ligand1,ligand2, '{2:},{3:}' poses1,poses2 '{4:}' str(MCSS).
        """
        
        init_command =  '$SCHRODINGER/run {0:}/main.py mcss INIT '.format(shared_paths['code'])
        init_command += '{0:} {1:} {2:} {3:} ' # ligand names
        init_command += self.init_file.format('{0:}-{1:}') + ' '
        init_command += self.atom_types + ' '
        init_command += '\n'
        
        rmsd_command = '$SCHRODINGER/run {0:}/main.py mcss RMSD '.format(shared_paths['code'])
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
        if l1 > l2: l1, l2 = l2, l1
        mcss = self.MCSSs["{}-{}".format(l1, l2)]
        if not mcss.is_valid(): return None
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
                if self.MCSSs[name].is_valid():
                    poseviewer_paths = {l1: self.pv_template.format(l1),
                                        l2: self.pv_template.format(l2)}
                    self.MCSSs[name].load_rmsds(self.rmsd_file.format(name),
                                                poseviewer_paths, max_poses)

    def verify_rmsds(self):
        """
        Check that all existing RMSD files are valid
        and contain entries for at least max_poses poses
        or the number of poses that exist for the given
        ligand. If ligands is None, check all pairs.

        * Delete all non-valid files *

        max_poses: int
        ligands: iterable
        """
        self.load_mcss()
        for name, mcss in self.MCSSs.items():
            rmsd_file_exists = os.path.exists(self.rmsd_file.format(name))
            if rmsd_file_exists and mcss.is_valid():
                poseviewer_paths = {mcss.l1: self.pv_template.format(mcss.l1),
                                    mcss.l2: self.pv_template.format(mcss.l2)}
                if not self.MCSSs[name].verify_rmsds(self.rmsd_file.format(name),
                                                     poseviewer_paths, max_poses):
                    print("Deleting RMSD file {}.".format(self.rmsd_file.format(name)))
                    os.system('rm {}'.format(self.rmsd_file.format(name)))
            self.MCSSs[name].rmsds = {} # To limit memory usage.

    def sort_by_mcss(self, query, ligands):
        """
        Sort ligands by size of maximum common substructure with query.
        Will give an error if any pairs have not been computed.

        * Must first call load_mcss *

        query: string, ligand name
        ligands: list(ligands), set of ligands to be sorted
        """
        self.load_mcss()
        def size(ligand):
            if ligand < query:
                name = '{}-{}'.format(ligand, query)
            else:
                name = '{}-{}'.format(query, ligand)
            return self.MCSSs[name].n_mcss_atoms
        return sorted(ligands, key=size, reverse=True)

    def get_mcss_size(self, l1, l2):
        """
        Returns the fraction of the smaller of the ligands
        that composes the MCSS for l1 and l2.

        l1, l2: string, ligand names
        """
        if l1 > l2: l1, l2 = l2, l1
        mcss = self.MCSSs["{}-{}".format(l1, l2)]
        return mcss.n_mcss_atoms / float(min(mcss.n_l1_atoms, mcss.n_l2_atoms))

    def load_mcss(self, temp_init_files = None):
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
            with open(temp) as fp:
                try:
                    mcss = MCSS.from_string(fp.readline().strip())
                except ValueError:
                    assert False, temp
                if mcss is not None:
                    self.MCSSs[mcss.name] = mcss

    # All of below are relevant for computation only
    def compute_mcss(self, chembl = True, pick_helpers={}):
        """
        Compute unfinished MCSS features. See above class description for more detail.

        chembl (bool): if True compute MCSS for PDB + CHEMBL ligand pairs,
                       else only consider PDB ligand pairs.
        pick_helpers: {pick_helpers_filename: {query_pdb: [chembl, ...]}}
        """
        previous_cwd = os.getcwd()
        os.system('mkdir -p {}'.format(self.root))
        os.chdir(self.root)

        print("{} PDB ligands, {} CHEMBL ligands".format(len(self.pdb),
                                                         len(self.chembl)))

        if pick_helpers:
            num_pdb = [len(v) for v in pick_helpers.values()]
            num_chembl = [len(v)
                          for _v in pick_helpers.values()
                          for v in _v.values()]
            if min(num_pdb) != max(num_pdb):
                print("# PDB ranges from {} to {}".format(min(num_pdb), max(num_pdb)))
            if min(num_chembl) != max(num_chembl):
                print("# CHEMBL ranges from {} to {}".format(min(num_chembl), max(num_chembl)))
            num_pdb = num_pdb[0]
            num_chembl = num_chembl[0]
            print("{} sort schemes in pick_helpers, "
                  "{} PDB ligands each, "
                  "{} CHEMBL ligands per PDB ligand".format(len(pick_helpers),
                                                            num_pdb, num_chembl))

        self._collate_mcss()
        self._add_pdb_to_pdb()
        self._add_crystal()
        if chembl:
            self._add_pdb_to_allchembl()
            self._add_pick_helpers(pick_helpers)
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
                self._add(l1, l2, compute_rmsds, compute_small = True)

    def _add_crystal(self):
        compute_rmsds = True
        crystal_lig = self.st + '_crystal_lig'
        for ligand in self.pdb:
            if ligand == self.st + '_lig': continue 
            self._add(crystal_lig, ligand, compute_rmsds)

    def _add_pdb_to_allchembl(self):
        """
        Add all pdb - chembl ligand pairs.
        Most of these are only used for deciding which chembl ligands to use
        so don't compute rmsd for all of them.
        """
        compute_rmsds = False
        for l1 in self.pdb:
            for l2 in self.chembl:
                self._add(l1,l2, compute_rmsds)

    def _add_pick_helpers(self, pick_helpers):
        """
        Adds pairs of pdb to chembl ligands and chembl to chembl ligands
        that will be jointly used in a score computation as specified
        by pick_helpers.

        Chembl - Chembl pairs that are specified by "chembl"
        pick_helpers: {pick_helpers_filename: {query_pdb: [chembl, ...]}}
        """
        compute_rmsds = True
        crystal_lig = self.st + '_crystal_lig'
        for fname, queries in pick_helpers.items():
            for query, chembl_ligands in queries.items():
                for i, l1 in enumerate(chembl_ligands):
                    assert l1 != '', pick_helpers # switch to continue if this happens
                    self._add(query,       l1, compute_rmsds)
                    self._add(crystal_lig, l1, compute_rmsds)
                    for l2 in chembl_ligands[i+1:]:
                        assert l2 != '', pick_helpers # switch to continue if this happens
                        self._add(l1, l2, compute_rmsds)

    def _add(self, l1, l2, compute_rmsd, compute_small = False):
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
              and (l1 in self.docked or 'crystal' in l1)
              and (l2 in self.docked or 'crsytal' in l2)
              and self.MCSSs[name].is_valid()
              and not os.path.exists(self.rmsd_file.format(name))):
            self.no_rmsd.add((l1, l2))
        elif (compute_small
              and not self.MCSSs[name].n_mcss_atoms
              and min(self.MCSSs[name].n_l1_atoms, self.MCSSs[name].n_l2_atoms) < 25
              and not self.MCSSs[name].tried_small):
            self.no_mcss_small.add((l1, l2))

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
        if self.no_mcss_small:
            print(len(self.no_mcss_small), 'small mcss init pairs left')
            self._execute_init(small = True)

    def _execute_init(self, small = False):
        for i, pairs in enumerate(grouper(self.INIT_GROUP_SIZE,
                                          self.no_mcss_small if small else self.no_mcss)):
            script = 'init{}.sh'.format(i)
            contents = 'export SCHRODINGER_CANVAS_MAX_MEM=1e+12\n'
            for l1, l2 in pairs:
                poses1 = self.lig_template.format(l1.replace('_crystal', ''))
                poses2 = self.lig_template.format(l2.replace('_crystal', ''))
                contents += self.init_command.format(l1, l2, poses1, poses2)
                if small:
                    contents = contents[:-1] # remove new line
                    contents += ' small\n'
            with open(script, 'w') as f:
                f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script))

    def _execute_rmsd(self):
        for i, pairs in enumerate(grouper(self.RMSD_GROUP_SIZE, self.no_rmsd)):
            script = 'rmsd{}.sh'.format(i)
            contents = ''
            for l1, l2 in pairs:
                mcss = self.MCSSs["{}-{}".format(l1, l2)]
                poses1 = (self.crystal_template.format(l1.replace('_crystal', ''))
                          if 'crystal' in l1 else
                          self.pv_template.format(l1))
                poses2 = (self.crystal_template.format(l2.replace('_crystal', ''))
                          if 'crystal' in l2 else
                          self.pv_template.format(l2))
                contents += self.rmsd_command.format(l1, l2, poses1, poses2, str(mcss))
            with open(script, 'w') as f:
                f.write(self.TEMPLATE.format(contents))
            os.system('sbatch {}'.format(script))

import os
from schrodinger.structure import StructureReader, StructureWriter

class MCSS:
    """
    Reads and writes MCSS features for a ligand pair.

    There are two key phases of the computation:
        (1) Identification of maximum common substructure(s)
        (2) Computation of RMSDs between the substructures in
            docking results.

    Task (1) is accomplished using Schrodinger's canvasMCSS utility.
    Task (2) is accomplished by identifying all matches of the substructure(s)
    from (1) and finding the pair with the mimimum RMSD. This is a subtlely
    difficult task because of symmetry concerns and details of extracting
    substructures.

    MCSS must be at least half the size of the smaller of the ligands
    or no RMSDs are computed.

    A key design decision is to not specify any file names in this class
    (other than those associated with temp files). The implication of this
    is that MCSSController will be completely in control of this task, while
    this class can be dedicated to actually computing the MCSS feature.
    """
    
    mcss_cmd = ("$SCHRODINGER/utilities/canvasMCS -imae {} -ocsv {}"
                " -stop {} -atomtype C {}")
    
    def __init__(self, l1, l2):
        """
        l1, l2: string, ligand names
        """
        if l1 > l2: l1, l2 = l2, l1

        self.l1 = l1
        self.l2 = l2
        self.name = "{}-{}".format(l1, l2)
        
        self.n_l1_atoms = 0
        self.n_l2_atoms = 0
        self.n_mcss_atoms = 0
        self.smarts_l1 = []
        self.smarts_l2 = []
        self.rmsds = {}

        self.tried_small = False

        # Deprecated.
        self.n_mcss_bonds = 0

    def __str__(self):
        return ','.join(map(str,
            [self.l1, self.l2,
            self.n_l1_atoms, self.n_l2_atoms, self.n_mcss_atoms, self.n_mcss_bonds,
            ';'.join(self.smarts_l1), ';'.join(self.smarts_l2), self.tried_small]
            ))

    def is_valid(self):
        """
        Decides if the mcss is large enough to include in calculations.
        
        (It would be problematic if small MCSSs were taken into account
        because the score is based solely on RMSD).
        """
        return (2 * self.n_mcss_atoms > min(self.n_l1_atoms, self.n_l2_atoms)
                and self.n_mcss_atoms > 10)
    
    # Constructors
    @classmethod
    def from_string(cls, S):
        """
        Creates an MCSS instance from a string returned by the __str__ method.
        """
        tok = S.strip().split(',')
        (l1, l2, n_l1_atoms, n_l2_atoms,
         n_mcss_atoms, n_mcss_bonds, smarts_l1, smarts_l2) = tok[:8]

        mcss = MCSS(l1, l2)
        mcss.n_l1_atoms = int(n_l1_atoms)
        mcss.n_l2_atoms = int(n_l2_atoms)
        mcss.n_mcss_atoms = int(n_mcss_atoms)
        mcss.n_mcss_bonds = int(n_mcss_bonds)
        mcss.smarts_l1 = smarts_l1.split(';')
        mcss.smarts_l2 = smarts_l2.split(';')

        tried_small = len(tok)==9 and tok[8] == 'True'
        mcss.tried_small = tried_small
        return mcss

    # MCSS Computation Methods.
    def compute_mcss(self, ligands, init_file, mcss_types_file, small = False):
        """
        Compute the MCSS file by calling Schrodinger canvasMCSS.

        Updates instance with MCSSs present in the file
        """
        structure_file = '{}.ligands.mae'.format(init_file)
        mcss_file = '{}.mcss.csv'.format(init_file)
        stwr = StructureWriter(structure_file)
        stwr.append(ligands[self.l1])
        stwr.append(ligands[self.l2])
        stwr.close()
        self._set_ligand_sizes(structure_file)
        
        if os.system(self.mcss_cmd.format(structure_file,
                                          mcss_file,
                                          5 if small else 10,
                                          mcss_types_file)):
            assert False, 'MCSS computation failed'
        self._set_mcss(mcss_file)
        self.tried_small = small

        with open(init_file, 'w') as fp:
            fp.write(str(self)+'\n')

        os.system('rm {} {}'.format(structure_file, mcss_file))

    def _set_ligand_sizes(self, structure_file):
        try:
            refs = [st for st in StructureReader(structure_file)]
        except:
            print('Unable to read MCSS structure file for', self.l1, self.l2)
            return None
        if len(refs) != 2:
            print('Wrong number of structures', self.l1, self.l2)
            return None
        ref1, ref2 = refs
        n_l1_atoms = len([a for a in ref1.atom if a.element != 'H'])
        n_l2_atoms = len([a for a in ref2.atom if a.element != 'H'])
        if self.n_l1_atoms:
            assert self.n_l1_atoms == n_l1_atoms
        if self.n_l2_atoms:
            assert self.n_l2_atoms == n_l2_atoms
        self.n_l1_atoms = n_l1_atoms
        self.n_l2_atoms = n_l2_atoms

    def _set_mcss(self, mcss_file):
        """
        Updates MCS from the direct output of canvasMCSS.

        Note that there can be multiple maximum common substructures
        of the same size.
        """
        ligs = {}
        n_mcss_atoms = None
        with open(mcss_file) as fp:
            fp.readline() # Header
            for line in fp:
                smiles, lig, _, _, _, _n_mcss_atoms, _n_mcss_bonds = line.strip().split(',')[:7]
                smarts = line.strip().split(',')[-1] # There are commas in some of the fields
                _n_mcss_atoms = int(_n_mcss_atoms)
                
                assert n_mcss_atoms is None or n_mcss_atoms == _n_mcss_atoms, self.name

                if lig not in ligs: ligs[lig] = []
                ligs[lig] += [smarts]
                n_mcss_atoms = _n_mcss_atoms

        if len(ligs) != 2:
            print('Wrong number of ligands in MCSS file', ligs)
            return None
        assert all(smarts for smarts in ligs.values()), ligs

        # MCSS size can change when tautomers change. One particularly prevalent
        # case is when oxyanions are neutralized. Oxyanions are sometimes specified
        # by the smiles string, but nevertheless glide neutralizes them.
        # Can consider initially considering oxyanions and ketones interchangable
        # (see mcss15.typ).
        if self.n_mcss_atoms:
            assert self.n_mcss_atoms <= n_mcss_atoms+1, 'MCSS size decreased by more than 1.'
            if self.n_mcss_atoms < n_mcss_atoms:
                print(self.name, 'MCSS size increased.')
            if self.n_mcss_atoms > n_mcss_atoms:
                print(self.name, 'MCSS size dencreased by one.')
        
        self.n_mcss_atoms = n_mcss_atoms
        self.smarts_l1 += ligs[self.l1.replace('_crystal', '')]
        self.smarts_l2 += ligs[self.l2.replace('_crystal', '')]

    # RMSD Methods.
    def load_rmsds(self, rmsd_file, poseviewer_paths, max_poses):
        """
        Loads RMSDs from rmsd_file for the top max_poses poses. 

        rmsd_file: string, absolute path to computed rmsds
        poseviewer_paths: {l1: poseviewer1, l2:poseviewer2},
            must contain self.l1 and self.l2
            paths to docking output for the two ligands
        max_poses: int, maximum number of poses to consider
        """
        with open(rmsd_file) as f:
            for line in f:
                line = line.strip().split(',')
                p1, p2 = int(line[0]), int(line[1])
                rmsd = float(line[2])
                if p1 < max_poses and p2 < max_poses:
                    self.rmsds[(p1,p2)] = rmsd

    def verify_rmsds(self, rmsd_file, poseviewer_paths, max_poses):
        """
        Verifies that all RMSDs have been computed (based on the number of
        poses in pose_viewer_paths).

        Returns False if RMSD file is not valid.
        """
        if not (os.path.exists(poseviewer_paths[self.l1])
            and os.path.exists(poseviewer_paths[self.l2])):
            return True
        if not os.path.exists(rmsd_file): return False
        
        from schrodinger.structure import StructureReader
        n_poses_l1 = len(list(StructureReader(poseviewer_paths[self.l1]))) - 1
        n_poses_l2 = len(list(StructureReader(poseviewer_paths[self.l2]))) - 1
        
        self.load_rmsds(rmsd_file, poseviewer_paths, max_poses)
        for i in range(min(max_poses, n_poses_l1)):
            for j in range(min(max_poses, n_poses_l2)):
                if (i, j) not in self.rmsds:
                    return False
        return True

    def write_rmsds(self, poseviewer_paths, init_file, mcss_types_file,
                    rmsd_file, max_poses): #, full_struct=False):
        """
        rmsd_file (str): path to where to write RMSDs.
        init_file (str): path to where to write MCSS init file if on-the-fly
            computation is necessary.
        pose_viewer_paths {self.l1: l1 docking, self.l2: l2 docking}
        max_poses (int): compute RMSDs for at most max_poses poses

        Computes the RMSD between MCSSs for all pairs of poses.
        """
        pv1 = list(StructureReader(poseviewer_paths[self.l1])
                   )[poseviewer_paths[self.l1][-8:] == 'pv.maegz':]
        pv2 = list(StructureReader(poseviewer_paths[self.l2])
                   )[poseviewer_paths[self.l2][-8:] == 'pv.maegz':]
        
        # These objects are lists of lists of lists of atom indices!
        l1_atom_idxss, l2_atom_idxss = self._get_atom_idxss(pv1[0], pv2[0],
                                                            init_file, mcss_types_file)
        
        # l1_atom_idxss = [atom.idx for atom in pv1[0].atom if atom.element != 'H']
        # l2_atom_idxss = [atom.idx for atom in pv2[0].atom if atom.element != 'H']

        rmsds = {}
        for i, pose1 in enumerate(pv1[:max_poses]):
            for j, pose2 in enumerate(pv2[:max_poses]):
                rmsd = float('inf')
                for l1_atom_idxs, l2_atom_idxs in zip(l1_atom_idxss, l2_atom_idxss):
                    for l1_atom_idx in l1_atom_idxs:
                        for l2_atom_idx in l2_atom_idxs:
                            rmsd = min(rmsd,
                                       self._calculate_rmsd(pose1, pose2,
                                                            l1_atom_idx, l2_atom_idx,
                                                            'mcss16' in mcss_types_file))
                if rmsd == float('inf'):
                    print("no mcss found"+','.join([str(i),str(j), self.l1,self.l2]))
                    print(l1_atom_idxss, l2_atom_idxss)
                    substructure1 = pose1.extract(l1_atom_idxss[0][0])
                    substructure2 = pose2.extract(l2_atom_idxss[0][0])
                    print(len(substructure1.atom), len(substructure2.atom))
                    print(len(substructure1.bond), len(substructure2.bond))
                    assert False
                rmsds[(i, j)] = rmsd

        with open(rmsd_file, 'w') as f:
            for (i, j), rmsd in rmsds.items():
                f.write('{},{},{}\n'.format(i, j, rmsd))

    def _get_atom_idxss(self, pose1, pose2, init_file, mcss_types_file):
        """
        Return atom indices of MCSS in pose1 and pose2.
        If initially, there is no match, try recomputing MCSS
        using docked pose.
        """

        # evaluate_smarts returns [[atom_index, ...], ...]
        l1_atom_idxss = [evaluate_smarts_canvas(pose1, smarts)
                         for smarts in self.smarts_l1]
        l2_atom_idxss = [evaluate_smarts_canvas(pose2, smarts)
                         for smarts in self.smarts_l2]

        # If initially no matches, recompute MCSSs using docked ligands.
        # This can occur do to changes in the tautomers produced by glide.
        if not any(len(l1_atom_idxs)*len(l2_atom_idxs)
                   for l1_atom_idxs, l2_atom_idxs in zip(l1_atom_idxss, l2_atom_idxss)):
            print('Recomputing mcss for {}-{}'.format(self.l1, self.l2))
            self.compute_mcss({self.l1: pose1, self.l2: pose2}, init_file, mcss_types_file)
            
            l1_atom_idxss = [evaluate_smarts_canvas(pose1, smarts)
                             for smarts in self.smarts_l1]
            l2_atom_idxss = [evaluate_smarts_canvas(pose2, smarts)
                             for smarts in self.smarts_l2]

        # If still no common substructure found, exit.
        assert len(l1_atom_idxss) and len(l2_atom_idxss), (l1_atom_idxss, l2_atom_idxss)
        assert any(len(l1_atom_idxs)*len(l2_atom_idxs)
                   for l1_atom_idxs, l2_atom_idxs in zip(l1_atom_idxss, l2_atom_idxss)), \
               (l1_atom_idxss, l2_atom_idxss, self.smarts_l1, self.smarts_l2)
        
        return l1_atom_idxss, l2_atom_idxss

    def _calculate_rmsd(self, pose1, pose2, atom_idx1, atom_idx2, merge_halogens):
        """
        Calculates the RMSD between the atoms atom_idx1 in pose1
        and the atoms atom_idx2 in pose2.
        
        pose1, pose2: schrodinger.structure
        atom_idx1, atom_idx2: [int, ...]
        merge_halogens: If true then change the atomic number of all halogens
                        to 9 (the atomic number of flourine) before computing
                        rmsds. This allows for MCSS that treat all halogens
                        the same.
        """
        substructure1 = pose1.extract(atom_idx1)
        substructure2 = pose2.extract(atom_idx2)
        if merge_halogens:
            self._merge_halogens(substructure1)
            self._merge_halogens(substructure2)
        try:
            calc = ConformerRmsd(substructure1, substructure2)
            calc.use_heavy_atom_graph = True
            rmsd = calc.calculate()
        except:
            # This is necessary because there is a bug in the
            # Schrodinger software that results in incorrect
            # atom indices being used when the heavy_atom_graph
            # is used. That being said, the above is more reliable
            # than the below, so should be tried first.
            calc = ConformerRmsd(substructure1, substructure2)
            calc.use_heavy_atom_graph = False
            rmsd = calc.calculate()
        return rmsd

    def _merge_halogens(self, structure):
        """
        Sets atomic number for all halogens to be that for flourine.
        This enable use of ConformerRmsd for atom typing schemes that
        merge halogens.
        """
        for atom in structure.atom:
            if atom.atomic_number in [9, 17, 35, 53]:
                atom.atomic_number = 9
            
def main(args):
    mode = args[1]
    if mode == 'INIT':
        l1, l2, l1_path, l2_path, init_file, mcss_types_file = args[2:8]
        small = len(args) == 9

        mcss = MCSS(l1, l2)
        with StructureReader(l1_path) as ligand1, StructureReader(l2_path) as ligand2:
            ligands = {l1: next(ligand1), l2: next(ligand2)}
            mcss.compute_mcss(ligands, init_file, mcss_types_file, small)
    elif mode == 'RMSD':
        from schrodinger.structutils.analyze import evaluate_smarts_canvas
        from schrodinger.structutils.rmsd import ConformerRmsd
        (l1, l2, pv1_path, pv2_path, init_file, mcss_types_file,
         rmsd_file, max_poses, mcss_string_rep) = args[2:]
        max_poses = int(max_poses)
        poseviewer_paths = {l1:pv1_path, l2:pv2_path}
        mcss = MCSS.from_string(mcss_string_rep)
        mcss.write_rmsds(poseviewer_paths, init_file, mcss_types_file,
                         rmsd_file, max_poses)
    else:
        assert False

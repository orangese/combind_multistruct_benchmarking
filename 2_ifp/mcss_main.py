import os
import sys


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

    Inside the code we stipulate that the MCSS be at least half the size of
    the smaller of the ligands, or no RMSDs are computed.

    A key design decision is to not specify any file names in this class
    (other than those associated with temp files). The implication of this
    is that MCSSController will be completely in control of this task, while
    this class can be dedicated to actually computing the MCSS feature.
    """
    
    structure_path = 'mcss_in.mae'
    canvasMCSS_csv = 'mcss.csv'
    mcss_cmd = ("$SCHRODINGER/utilities/canvasMCS -imae {} -ocsv {}"
                " -stop 10 -atomtype C {}").format(structure_path, canvasMCSS_csv, '{}')
    
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
        self.n_mcss_bonds = 0
        self.smarts_l1 = []
        self.smarts_l2 = []
        
        self.rmsds = {}

        # Cache for atom indexes of MCSS matches in pose viewers
        self.atom_idx_l1 = None
        self.atom_idx_l2 = None

    def __str__(self):
        return ','.join(map(str,
            [self.l1, self.l2,
            self.n_l1_atoms, self.n_l2_atoms, self.n_mcss_atoms, self.n_mcss_bonds,
            ';'.join(self.smarts_l1), ';'.join(self.smarts_l2)]
            ))

    def is_valid(self):
        """
        Decides if the mcss is large enough to include in calculations.
        (It would be problematic if small MCSSs were taken into account
        because the score is based solely on RMSD).

        * Important! *
        """
        return 2 * self.n_mcss_atoms > min(self.n_l1_atoms, self.n_l2_atoms)
    
    # Constructors
    @classmethod
    def from_string(cls, S):
        """
        Creates an MCSS instance from a string returned by the __str__ method.
        """
        (l1, l2, n_l1_atoms, n_l2_atoms,
         n_mcss_atoms, n_mcss_bonds, smarts_l1, smarts_l2) = S.strip().split(',')

        mcss = MCSS(l1, l2)
        mcss.n_l1_atoms = int(n_l1_atoms)
        mcss.n_l2_atoms = int(n_l2_atoms)
        mcss.n_mcss_atoms = int(n_mcss_atoms)
        mcss.n_mcss_bonds = int(n_mcss_bonds)
        mcss.smarts_l1 = smarts_l1.split(';')
        mcss.smarts_l2 = smarts_l2.split(';')
        return mcss

    @classmethod
    def from_canvasMCSS(cls, dir_name):
        """
        Reads the direct output of canvasMCSS and returns an MCSS instance.
        """
        from schrodinger.structure import StructureReader
        mcss_output = "{}/{}".format(dir_name, cls.canvasMCSS_csv)
        mcss_struct = "{}/{}".format(dir_name, cls.structure_path)
        if not (os.path.exists(mcss_output) and os.path.exists(mcss_struct)): return None
        try:
            refs = [st for st in StructureReader(mcss_struct)]
        except:
            print('Unable to read MCSS structure file for', dir_name)
            return None
        if len(refs) != 2:
            print('Wrong number of structures', dir_name)
            return None
        ligs = {}
        n_mcss_bonds = None
        n_mcss_atoms = None
        with open(mcss_output) as f:
            f.readline() # Header
            for line in f:
                smiles, lig, _, _, _, _n_mcss_atoms, _n_mcss_bonds = line.strip().split(',')[:7]
                smarts = line.strip().split(',')[-1] # There are commas in some of the fields
                _n_mcss_atoms, _n_mcss_bons = int(_n_mcss_atoms), int(_n_mcss_bonds)
                
                assert n_mcss_atoms is None or n_mcss_atoms == _n_mcss_atoms, dir_name
                assert n_mcss_bonds is None or n_mcss_bonds == _n_mcss_bonds, dir_name
                assert lig in dir_name, dir_name

                if lig not in ligs: ligs[lig] = []
                ligs[lig] += [smarts]
                n_mcss_atoms = _n_mcss_atoms
                n_mcss_bonds = _n_mcss_bonds
        if len(ligs) != 2:
            print('Wrong number of ligands in MCSS file', ligs, dir_name)
            return None
        #assert len(ligs) == 2, dir_name
        assert all(smarts for smarts in ligs.values()), dir_name

        mcss = MCSS(*ligs.keys())
        mcss.n_mcss_atoms = n_mcss_atoms
        mcss.n_mcss_bonds = n_mcss_bonds
        mcss.smarts_l1 = ligs[mcss.l1]
        mcss.smarts_l2 = ligs[mcss.l2]

        # Get ligand sizes.
        # The below line can fail if there are multiple compies of each ligand
        # writtent of the mcss_struct file
        ref1, ref2 = refs
        mcss.n_l1_atoms = len([a for a in ref1.atom if a.element != 'H'])
        mcss.n_l2_atoms = len([a for a in ref2.atom if a.element != 'H'])
        return mcss

    def compute_mcss(self, mcss_types_file, ligand_paths):
        """
        Compute the MCSS file by calling Schrodinger canvasMCSS.

        Copies ligands to self.structure_path with self.l1 first
        in the file.

        Writes MCSS to self.canvasMCSS_csv.

        * Must be called from MCSS temp directory *
        """
        if os.path.exists(self.structure_path):
            os.system("rm {}".format(self.structure_path))
        stwr = StructureWriter(self.structure_path)
        stwr.append(next(StructureReader(ligand_paths[self.l1])))
        stwr.append(next(StructureReader(ligand_paths[self.l2])))
        stwr.close()
        
        os.system(self.mcss_cmd.format(mcss_types_file))

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
        self._verify_rmsds(rmsd_file, poseviewer_paths, max_poses)

    def _verify_rmsds(self, rmsd_file, poseviewer_paths, max_poses):
        """
        Verifies that all RMSDs have been computed (based on the number of
        poses in pose_viewer_paths).

        Returns False if RMSD file is not valid.
        """
        if not (os.path.exists(poseviewer_paths[self.l1])
            and os.path.exists(poseviewer_paths[self.l2])):
            return True
        if not os.path.exists(rmsd_file): return False

        n_poses_l1 = len(StructureReader(poseviewer_paths[self.l1])) - 1
        n_poses_l2 = len(StructureReader(poseviewer_paths[self.l2])) - 1
        for i in min(max_poses, n_poses_l1):
            for j in min(max_poses, n_poses_l2):
                assert (i, j) in self.rmsds

    def write_rmsds(self, rmsd_file, poseviewer_paths, max_poses):
        """
        Computes the RMSD between MCSSs for all pairs of poses.
        """
        pv1 = list(StructureReader(poseviewer_paths[self.l1]))[1:]
        pv2 = list(StructureReader(poseviewer_paths[self.l2]))[1:]
        # evaluate_smarts returns [[atom_index, ...], ...]
        self.atom_idx_l1 = [evaluate_smarts(pv1[0], smarts, unique_sets=True)
                            for smarts in self.smarts_l1]
        self.atom_idx_l2 = [evaluate_smarts(pv2[0], smarts, unique_sets=True)
                            for smarts in self.smarts_l2]

        pv1 = [self._unproc_st(pose) for pose in pv1]
        pv2 = [self._unproc_st(pose) for pose in pv2]
        with open(rmsd_file, 'w') as f:
            for i, p1 in enumerate(pv1[:max_poses]):
                for j, p2 in enumerate(pv2[:max_poses]):
                    rmsd = float('inf')
                    for m1, m2 in self._find_mcss_matches(p1, p2):
                        conf_rmsd = ConformerRmsd(m1, m2)
                        conf_rmsd.use_heavy_atom_graph = True
                        rmsd = min(rmsd, conf_rmsd.calculate())

                    assert rmsd != float('inf'), "no mcss found"+','.join(str(i),str(j),
                        self.l1,self.l2)
                    f.write('{},{},{}\n'.format(i, j, rmsd))

    # Utilities for matching atoms in MCSS
    def _find_mcss_matches(self, st1, st2):
        """
        Finds all matches between pairs of substructures
        specificed by self.smarts in st1 and st2.
        Returns list of tuples of schrodinger.structure.Structure
        objects where the first entry is for st1 and second is for st2.

        st1, st2: schrodinger.structure.Structure, ligands
        """
        all_pairs = []
        for mcss1 in self.atom_idx_l1:
            for mcss2 in self.atom_idx_l2:
                # st.extract returns a schrodinger.structure.Structure
                list1 = [st1.extract(m) for m in mcss1]
                list2 = [st2.extract(m) for m in mcss2]

                for m1 in list1:
                    for m2 in list2:
                        attempt_match = self._match(m1, m2)
                        if attempt_match is not None:
                            all_pairs.append((m1,m2))
        return all_pairs

    def _match(self, st1, st2):
        """
        Attempts to find a match between the two given structures.
        If no match immediately present, deletes extraneous bonds.
        Returns schrodinger.structure.Structure tuple of match if
        found otherwise None.

        st1, st2: schrodinger.structure.Structure, ligand substructures
        """
        if st1.isEquivalent(st2, False): return (st1, st2)
        if abs(len(st1.bond) - len(st2.bond)) == 1:
            return self._delete_extraneous_bonds(st1, st2)
        if len(st1.bond) == len(st2.bond) == self.n_mcss_bonds + 1:
            return self._delete_bond_pair(st1, st2)

    def _unproc_st(self, st):
        """
        Set bond orders to 1 and partial charges to 0.

        st1: schrodinger.structure.Structure
        """
        st2 = st.copy()
        for b in st2.bond: b.order = 1
        for a in st2.atom: a.formal_charge = 0
        st2.retype()
        return st2

    def _delete_extraneous_bonds(self, st1, st2):
        if len(st1.bond) + 1 == len(st2.bond): # st2 has the extra bond
            st1, st2 = st2, st1
        for bond in st1.bond:
            atom1,atom2 = bond.atom1, bond.atom2
            st1.deleteBond(atom1, atom2)
            if st1.isEquivalent(st2,False):
                return (st1,st2)
            st1.addBond(atom1, atom2, bond.order)

    def _delete_bond_pair(self, st1, st2):
        for bond1 in st1.bond:
            order1 = bond1.order
            atom11, atom12 = bond1.atom1, bond1.atom2
            st1.deleteBond(atom11,atom12)
            for bond2 in st2.bond:
                order2 = bond2.order
                atom21, atom22 = bond2.atom1, bond2.atom2
                st2.deleteBond(atom21, atom22)
                if st1.isEquivalent(st2,False):
                    return (st1,st2)
                st2.addBond(atom21, atom22, order2)
            st1.addBond(atom11, atom12, order1)

if __name__ == '__main__':
    from schrodinger.structure import StructureReader, StructureWriter
    mode = sys.argv[1]
    if mode == 'INIT':
        l1, l2, l1_path, l2_path, mcss_types_file = sys.argv[2:]
        mcss = MCSS(l1, l2)
        ligand_paths = {l1:l1_path, l2: l2_path}
        mcss.compute_mcss(mcss_types_file, {l1: l1_path, l2: l2_path})
    elif mode == 'RMSD':
        from schrodinger.structutils.analyze import evaluate_smarts
        from schrodinger.structutils.rmsd import ConformerRmsd
        l1, l2, pv1_path, pv2_path, rmsd_file, max_poses, mcss_string_rep = sys.argv[2:]
        max_poses = int(max_poses)
        poseviewer_paths = {l1:pv1_path, l2:pv2_path}
        mcss = MCSS.from_string(mcss_string_rep)
        mcss.write_rmsds(rmsd_file, poseviewer_paths, max_poses)
    else:
        assert False

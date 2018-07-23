import os
import sys

class MCSS:
    """
    Reads and writes MCSS features for a ligand pair.
    """
    def __init__(self, l1, l2, max_poses,
                 st, prot, mcss_version, docking_version, data_path, code_path):
        if l1 > l2: l1, l2 = l2, l1
        self.l1 = l1
        self.l2 = l2
        self.max_poses = max_poses

        self.mcss_size = -1
        self.ligand_size = {l1:0, l2:0}
        self.n_matches = 0
        self.smarts = {l1:[], l2:[]}
        self.n_bonds = -1
        
        self.rmsds = {}

        # Paths
        self.name = "{}-{}".format(self.l1, self.l2)
        self.root = "{}/{}".format(data_path, prot)
        
        # Processed ligands
        self.l1_ref_path = "{0:}/ligands/prepared_ligands/{1:}/{1:}.mae".format(self.root, self.l1)
        self.l2_ref_path = "{0:}/ligands/prepared_ligands/{1:}/{1:}.mae".format(self.root, self.l2)
        
        # Docking results
        self.l1_pv_path  = "{0:}/docking/{1:}/{2:}-to-{3:}/{2:}-to-{3:}_pv.maegz".format(
            self.root, docking_version, self.l1, st)
        self.l2_pv_path  = "{0:}/docking/{1:}/{2:}-to-{3:}/{2:}-to-{3:}_pv.maegz".format(
            self.root, docking_version, self.l2, st)
        
        # MCSS related files
        self.mcss_root = "{}/mcss/{}/{}".format(self.root, mcss_version, self.name)
        self.init_structure_path = 'mcss_in.mae'
        self.init_path = "{}.csv".format(self.name)
        self.size_path = "{}.size".format(self.name)
        self.rmsd_path = "{}-{}-{}.csv".format(self.name, st, docking_version)

        # Command to compute MCSS!
        mcss_type_file = "{}/2_ifp/custom_types/{}.typ".format(code_path, mcss_version)
        self.init_command = ("$SCHRODINGER/utilities/canvasMCS "
                             "-imae {} -ocsv {} -stop 10 "
                             "-atomtype C {}"
                             ).format(self.init_structure_path, self.init_path, mcss_type_file)

    def load_rmsds(self):
        """
        Load the RMSDs for this ligand pair. This is the primary method that should
        be run for accessing the MCSS features.
        """
        self._load_size_file()
        if self.n_matches > 0:
            self._load_rmsd_file()

    def get_size(self):
        """
        Returns: ([int, int], float), the ligand sizes and MCSS sizes.
        """
        self._load_size_file()
        return self.ligand_size.values(), self.mcss_size

    # 1. INIT
    def write_init_file(self):
        """
        Compute the MCSS file by calling Schrodinger canvasMCSS.
        """
        if os.path.exists(self.init_structure_path):
            os.system("rm {}".format(self.init_structure_path))
        stwr = StructureWriter(self.init_structure_path)
        stwr.append(StructureReader(self.l1_ref_path).next())
        stwr.append(StructureReader(self.l2_ref_path).next())
        stwr.close()

        os.system(self.init_command)
        os.system("rm {}".format(self.init_structure_path))

    def _load_init_file(self):
        """
        Reads the direct output of canvasMCSS and loads relevant
        into instance.

        Verifies that all of the data is properly read. Will error
        if file is unreadable.

        Sets
         - self.mcss_size
         - self.n_bonds
         - self.smarts
        """
        with open("{}/{}".format(self.mcss_root, self.init_path)) as f:
            f.readline() # Header
            for line in f:
                smiles, lig, _, _, _, mcss_size, n_bonds = line.strip().split(',')[:7]
                smarts = line.strip().split(',')[-1] # There are commas in some of the fields
                self.smarts[lig] += [smarts]
                assert self.mcss_size == -1 or self.mcss_size == int(mcss_size)
                assert self.n_bonds == -1 or self.n_bonds == int(n_bonds)
                self.mcss_size = int(mcss_size)
                self.n_bonds = int(n_bonds)
        assert self.l1 in self.smarts
        assert self.l2 in self.smarts

    # 2. SIZE
    def write_size_file(self):
        """
        Computes the size of each ligand and verifies that we can
        identify the MCSS.

        Decides on a minimum size of MCSS to consider.
        """
        self._load_init_file()
        ref1 = StructureReader(self.l1_ref_path).next()
        ref2 = StructureReader(self.l2_ref_path).next()
        
        # Get ligand sizes
        s1 = len([a for a in ref1.atom if a.element != 'H'])
        s2 = len([a for a in ref2.atom if a.element != 'H'])
        
        # Don't use MCSS if the common substructure is less than half the size of the ligand
        if self.mcss_size * 2 < min(s1, s2):
            n_matches = -1 
        else:
            n_matches = len(self._find_mcss_matches(ref1, ref2))
        
        with open("{}/{}".format(self.mcss_root, self.size_path),'w') as f:
            f.write('{} matches\n'.format(n_matches))
            for smarts in self.smarts[self.l1]:
                f.write(','.join([self.l1, str(s1), str(self.mcss_size), smarts])+'\n')
            for smarts in self.smarts[self.l2]:
                f.write(','.join([self.l2, str(s2), str(self.mcss_size), smarts])+'\n')

    def _load_size_file(self):
        with open("{}/{}".format(self.mcss_root, self.size_path)) as f:
            self.n_matches = int(f.readline().split(' ')[0])
            for i,line in enumerate(f):
                lig, lsize, msize, smarts = line.strip().split(',')
                self.mcss_size = int(msize)
                self.ligand_size[lig] = int(lsize)
                self.smarts[lig].append(smarts)

        # Verify
        assert self.n_matches == -1 or self.n_matches > 0, self.n_matches
        assert self.l1 in self.ligand_size
        assert self.l2 in self.ligand_size
        assert self.l1 in self.smarts
        assert self.l2 in self.smarts

    # 3. RMSD
    def write_rmsd_file(self):
        """
        Computes the RMSD between MCSSs for all pairs of poses.
        """
        self._load_size_file()
        pv1 = list(StructureReader(self.l1_pv_path))[1:]
        pv2 = list(StructureReader(self.l2_pv_path))[1:]
        with open("{}/{}".format(self.mcss_root, self.rmsd_path), 'w') as f:
            if self.n_matches < 1: return # Leave empty file to mark as tried
            for i, p1 in enumerate(pv1):
                if i > self.max_poses: continue
                for j, p2 in enumerate(pv2):
                    if j > self.max_poses: continue
                    rmsd = float('inf')
                    for m1, m2 in self._find_mcss_matches(p1,p2):
                        conf_rmsd = ConformerRmsd(m1, m2)
                        conf_rmsd.use_heavy_atom_graph = True
                        rmsd = min(rmsd, conf_rmsd.calculate())

                    assert rmsd != float('inf'), "no mcss found"+','.join(str(i),str(j),
                        self.l1,self.l2)
                    f.write('{},{},{}\n'.format(i, j, rmsd))

    def _load_rmsd_file(self):
        with open("{}/{}".format(self.mcss_root, self.rmsd_path)) as f:
            for i, line in enumerate(f):
                line = line.strip().split(',')
                p1, p2 = int(line[0]), int(line[1])
                rmsd = float(line[2])
                self.rmsds[(p1,p2)] = rmsd

        # Assert that all RMSDs have been computed
        n_poses1 = len(list(StructureReader(self.l1_pv_path))[1:])
        n_poses2 = len(list(StructureReader(self.l2_pv_path))[1:])
        for i in range(min(n_poses1, self.max_poses)):
            for j in range(min(n_poses1, self.max_poses)):
                assert (i, j) in self.rmsds, "{} {} missing".format(i, j)

    # Utilities for matching atoms in MCSS
    def _find_mcss_matches(self, st1, st2, unproc=True):
        """
        Finds all matches between pairs of substructures
        specificed by self.smarts in st1 and st2.
        Returns list of tuples of schrodinger.structure.Structure
        objects where the first entry is for st1 and second is for st2.

        st1, st2: schrodinger.structure.Structure, ligands
        """
        all_pairs = []
        for smarts1 in self.smarts[self.l1]:
            for smarts2 in self.smarts[self.l2]:
                # evaluate_smarts returns [[atom_index, ...], ...]
                mcss1 = evaluate_smarts(st1, smarts1, unique_sets=True)
                mcss2 = evaluate_smarts(st2, smarts2, unique_sets=True)
                # st.extract returns a schrodinger.structure.Structure
                if unproc:
                    list1 = [self._unproc_st(st1.extract(m)) for m in mcss1]
                    list2 = [self._unproc_st(st2.extract(m)) for m in mcss2]
                else:
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
        if len(st1.bond) == len(st2.bond) == self.n_bonds + 1:
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

    # The below probably occur in weird cases, going to leave
    # blocked off until these come up and I identify them.
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
    """
    All of these must be called from the ligand pair's
    subdirectory! (Not true for loading methods)
    """
    from schrodinger.structure import StructureReader, StructureWriter
    from schrodinger.structutils.analyze import evaluate_smarts
    from schrodinger.structutils.rmsd import ConformerRmsd
    mode, l1, l2, max_poses, st, prot, mcss, docking, data, code = sys.argv[1:]
    mcss = MCSS(l1, l2, int(max_poses), st, prot,  mcss, docking, data, code)
    if mode == 'INIT':
        mcss.write_init_file()
    elif mode == 'SIZE':
        mcss.write_size_file()
    elif mode == 'RMSD':
        mcss.write_rmsd_file()
    else:
        assert False, ("Unclear what the purpose of this code is,"
                       "if turns out to be re-add from prior commit")

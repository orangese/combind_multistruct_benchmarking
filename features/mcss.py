import tempfile
import numpy as np
import subprocess
import os
from schrodinger.structure import StructureReader, StructureWriter, SmilesStructure
from schrodinger.structutils.rmsd import ConformerRmsd
from schrodinger.structutils.analyze import generate_smiles
from schrodinger.structutils.analyze import evaluate_smarts_canvas

def mcss(pv1, pv2, mcss_types_file, max_poses):
    with StructureReader(pv1) as sts1, StructureReader(pv2) as sts2:
        next(sts1), next(sts2) # Skip receptor
        sts2 = [merge_halogens(st) for st in sts2][:max_poses]

        rmsds = []
        for j, st1 in enumerate(sts1):
            if j == max_poses:
                break
            st1 = merge_halogens(st1)
            rmsds += [np.zeros(len(sts2)) + float('inf')]

            for i, st2 in enumerate(sts2):
                if i == max_poses:
                    break
                mcss = compute_mcss(st1, st2, mcss_types_file)
                mcss_st = SmilesStructure(mcss['st1'][0].split(',')[0].upper()).get2dStructure()
                
                mcss_atoms = n_atoms(mcss_st)
                st1_atoms = n_atoms(st1)
                st2_atoms = n_atoms(st2)

                if (2*mcss_atoms <= min(st1_atoms, st2_atoms)
                    or mcss_atoms <= 10):
                    continue

                for smarts1, smarts2 in zip(mcss['st1'], mcss['st2']):
                    atom_idxs1 = evaluate_smarts_canvas(st1, smarts1)
                    atom_idxs2 = evaluate_smarts_canvas(st2, smarts2)

                    for atom_idx1 in atom_idxs1:
                        for atom_idx2 in atom_idxs2:
                            rmsd = calculate_rmsd(st1, st2, atom_idx1, atom_idx2)
                            rmsds[-1][i] = min(rmsds[-1][i], rmsd)
    return np.vstack(rmsds)

def compute_mcss(st1, st2, mcss_types_file, memo={}):
    smi1 = generate_smiles(st1)
    smi2 = generate_smiles(st2)
    if (smi1, smi2) not in memo:
        cmd = "$SCHRODINGER/utilities/canvasMCS -imae {} -ocsv {} -stop 10 -atomtype C {}"
        with tempfile.TemporaryDirectory() as wd:
            mae = wd+'/temp.maegz'
            csv = wd+'/temp.csv'

            st1.title = 'st1'
            st2.title = 'st2'
            stwr = StructureWriter(mae)
            stwr.append(st1)
            stwr.append(st2)
            stwr.close()

            r = subprocess.run(cmd.format(os.path.basename(mae), os.path.basename(csv),
                                          os.path.abspath(mcss_types_file)),
                               cwd=wd, shell=True, stderr=subprocess.PIPE)


            # mcss can fail with memory usage error, generally when macrocycles
            # are present. just skip such cases.
            failed = 'memory usage' in str(r.stderr)

            if not failed:
                assert os.path.exists(csv)
                memo[(smi1, smi2)] = {'st1': [], 'st2': []}
                with open(csv) as fp:
                    fp.readline()
                    for line in fp:
                        lig = line.strip().split(',')[1]
                        smarts = line.strip().split(',')[-1]
                        memo[(smi1, smi2)][lig] += [smarts]
            else:
                print('mcss failed', smi1, smi2)
                memo[(smi1, smi2)] = {'st1': ['C'], 'st2': ['C']}
    return memo[(smi1, smi2)]

def calculate_rmsd(pose1, pose2, atom_idx1, atom_idx2):
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
    try:
        calc = ConformerRmsd(substructure1, substructure2)
        calc.use_heavy_atom_graph = True
        calc.renumber_structures = True
        rmsd = calc.calculate()
    except:
        # This is necessary because there is a bug in the
        # Schrodinger software that results in incorrect
        # atom indices being used when the heavy_atom_graph
        # is used. That being said, the above is more reliable
        # than the below, so should be tried first.
        calc = ConformerRmsd(substructure1, substructure2)
        calc.use_heavy_atom_graph = True
        calc.renumber_structures = False
        rmsd = calc.calculate()
    return rmsd

def merge_halogens(structure):
    """
    Sets atomic number for all halogens to be that for flourine.
    This enable use of ConformerRmsd for atom typing schemes that
    merge halogens.
    """
    for atom in structure.atom:
        if atom.atomic_number in [9, 17, 35, 53]:
            atom.atomic_number = 9
    return structure

def n_atoms(st):
    return sum(atom.element != 'H' for atom in st.atom)

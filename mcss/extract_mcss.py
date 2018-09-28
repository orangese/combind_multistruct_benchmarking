from containers import LigandManager, Dataset
from shared_paths import shared_paths
import os
import sys
n_ligands = int(sys.argv[1])
n_poses = 100

os.chdir(shared_paths['data'])
with open('/scratch/PI/rondror/combind/bpp_outputs/mcss.csv', 'w') as out:
    for protein in [d for d in sorted(os.listdir('.')) if d[0] != '.']:
        if protein == 'ERA': continue
        os.chdir(protein)
        print(protein)

        lm = LigandManager(shared_paths, protein)
        ligands = lm.docked(lm.pdb)[:n_ligands+1]
        self_docked = lm.st+'_lig'
        if self_docked in ligands:
            ligands.remove(self_docked)
        else:
            ligands.pop(-1)

        dataset = Dataset(shared_paths, [protein])
        dataset.load({protein: ligands}, load_fp=False, load_mcss=True)

        lm = dataset.proteins[protein].lm
        
        for i, ligand1 in enumerate(ligands):
            poses1 = dataset.proteins[protein].docking[lm.st].ligands[ligand1].poses
            for ligand2 in ligands[i+1:]:
                poses2 = dataset.proteins[protein].docking[lm.st].ligands[ligand2].poses
                if lm.mcss.get_rmsd(ligand1, ligand2, 0, 0) is None:
                    out.write(','.join([protein,ligand1,ligand2,'None'])+'\n')
                for r1 in range(n_poses):
                    if lm.mcss.get_rmsd(ligand1, ligand2, r1, 0) is None:
                        continue
                    for r2 in range(n_poses):
                        rmsd = lm.mcss.get_rmsd(ligand1, ligand2, r1, r2)
                        if rmsd is None: continue
                        out.write(','.join(map(str, [protein,ligand1,ligand2, r1, r2, rmsd,
                                                     poses1[r1].rmsd, poses2[r2].rmsd]))+'\n')
        os.chdir(shared_paths['data'])

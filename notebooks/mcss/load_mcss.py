import sys
import os
sys.path.append('../../1_dock')
sys.path.append('../../3_analyze')
sys.path.append('../../mcss')
sys.path.append('../../ifp')
from containers import Protein
from shared_paths import shared_paths

def add_ligpair_mcss(lm, protein, ligand1, ligand2,
                     ligand_size1, ligand_size2, mcss_size, ligs,
                     rmsd1, rmsd2, gscore1, gscore2, mcss):
    m = lm.mcss.MCSSs['{}-{}'.format(ligand1.ligand, ligand2.ligand)]
    if not m.rmsds: return
    # Ligand Properties
    ligs += [(protein, ligand1.ligand, ligand2.ligand)]
    ligand_size1 += [m.n_l1_atoms]
    ligand_size2 += [m.n_l2_atoms]
    mcss_size += [m.n_mcss_atoms]

    # Pose Properties
    rmsd1 += [[]]
    rmsd2 += [[]]
    gscore1 += [[]]
    gscore2 += [[]]
    mcss += [[]]

    for r1, pose1 in enumerate(ligand1.poses[:100]):
        for r2, pose2 in enumerate(ligand2.poses[:100]):
            rmsd1[-1] += [pose1.rmsd]
            rmsd2[-1] += [pose2.rmsd]
            gscore1[-1] += [pose1.gscore]
            gscore2[-1] += [pose2.gscore]
            mcss[-1] += [m.rmsds[(r1, r2)]]

def load_mcss(crystal = False, max_proteins = 100, max_ligands = 20):
    # Pose Properties
    rmsd1, rmsd2, gscore1, gscore2, mcss = [], [], [], [], []
    
    # Ligand Properties
    ligand_size1, ligand_size2, mcss_size, ligs = [], [], [], []
    
    proteins = [fname for fname in os.listdir(shared_paths['data']) if fname[0] != '.']
    
    for protein in proteins[:max_proteins]:
        print protein
        prot = Protein(protein)
        lm = prot.lm
        ligands = lm.docked(lm.pdb)[:max_ligands+1]
        self_docked = lm.st+'_lig'
        if self_docked in ligands:
            ligands.remove(self_docked)
        else:
            ligands.pop(-1)

        if crystal:
            prot.load_docking(ligands)
            crystal_lig = lm.st + '_crystal_lig'
            prot.load_docking([crystal_lig], load_mcss = True, load_crystal = True)
            for name in ligands:
                ligand1 = prot.docking[lm.st].ligands[crystal_lig]
                ligand2 = prot.docking[lm.st].ligands[name]
                add_ligpair_mcss(lm, protein, ligand1, ligand2,
                                 ligand_size1, ligand_size2, mcss_size, ligs,
                                 rmsd1, rmsd2, gscore1, gscore2, mcss)
                
        else:
            prot.load_docking(ligands, load_mcss = True)

            for i, name1 in enumerate(ligands[:max_ligands]):
                for name2 in ligands[i+1:max_ligands]:
                    if name1 > name2: name1, name2 = name2, name1
                    ligand1 = prot.docking[lm.st].ligands[name1]
                    ligand2 = prot.docking[lm.st].ligands[name2]
                    add_ligpair_mcss(lm, protein, ligand1, ligand2,
                                     ligand_size1, ligand_size2, mcss_size, ligs,
                                     rmsd1, rmsd2, gscore1, gscore2, mcss)
    return (ligand_size1, ligand_size2, mcss_size, ligs,
            rmsd1, rmsd2, gscore1, gscore2, mcss)

def flatten(ligand_size1, ligand_size2, mcss_size, ligs,
            rmsd1, rmsd2, gscore1, gscore2, mcss):
    X = []
    for ls1, ls2, ms, _, r1, r2, g1, g2, m in zip(ligand_size1, ligand_size2, mcss_size, ligs,
                                                  rmsd1, rmsd2, gscore1, gscore2, mcss):
        X += [np.vstack([[ls1]*len(r1), [ls2]*len(r1), [ms]*len(r1),
                         r1, r2, g1, g2, m]).T]
    return np.vstack(X)

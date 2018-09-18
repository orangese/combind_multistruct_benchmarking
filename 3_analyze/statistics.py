import sys
import numpy as np
from density_estimate import DensityEstimate
from pairs import LigPair
from containers import Dataset

sys.path.append('../1_dock/')
from shared_paths import shared_paths

def merge_stats(stats1, stats2):
    '''
    stats {'native' or 'reference': {interaction (str): DensityEstimate}}
    
    Merges the dictionaries of DensityEstimate's stats1 and stats2.
    '''
    out = {}
    for d in set(list(stats1.keys()) + list(stats2.keys())):
        out[d] = {}
        I = []
        if d in stats1: I += list(stats1[d].keys())
        if d in stats2: I += list(stats2[d].keys())
        for i in set(I):
            if   d not in stats1 or i not in stats1[d]:
                out[d][i] = stats2[d][i]
            elif d not in stats2 or i not in stats2[d]:
                out[d][i] = stats1[d][i]
            else:
                out[d][i] = stats1[d][i].average(stats2[d][i])
    return out

def statistics_lig_pair(lig_pair, interactions, p_native, weighted, sd):
    stats = {'native':{}, 'reference':{}}
    for interaction in interactions:
        X_native, X_ref, w_ref = [], [], []
        for (r1,r2), pp in lig_pair.pose_pairs.items():
            pp_x = lig_pair.get_feature(interaction, r1, r2)
            if pp_x is not None:
                if pp.correct():
                    X_native += [pp_x]
                X_ref += [pp_x]
                w_ref += [lig_pair.get_gscores(r1, r2)]

        X_native, X_ref = np.array(X_native), np.array(X_ref)
        w_ref = np.array([p_native(g1) * p_native(g2) for (g1, g2) in w_ref])
        if not weighted: w_ref = 1

        stats['native'][interaction]    = DensityEstimate(domain = (0, 1),
                                                          sd = sd,
                                                          reflect = True)
        stats['native'][interaction].fit(X_native)
        stats['reference'][interaction] = DensityEstimate(domain = (0, 1),
                                                          sd = sd,
                                                          reflect = True)
        stats['reference'][interaction].fit(X_ref, w_ref)
    return stats

def statistics(data, interactions, pnative, points=100,
               weighted = True, mcss=True, max_poses = 100, sd = 0.005):
    '''
    data    {protein (str): [ligand (str), ...]}
    interactions {name (str): [code (int), ...]}
    pnative (DensityEstimate): Likelihood of a pose being correct
        given its glide score.

    Returns DensityEstimate's representing the native and reference
    distribution of the statistics for all proteins, ligands in data.

    For each ligand pair, write distribution to file after computed and
    subsequently just read from this file.
    '''
    stats = {}
    for protein, ligands in data.items():
        print(protein, len(ligands))
        protein_stats = {}
        dataset = Dataset(shared_paths, [protein])
        dataset.load({protein: ligands}, load_fp=True, load_mcss=mcss)
        st = dataset.proteins[protein].lm.st
        docking = dataset.proteins[protein].docking[st]
        for i, ligand1 in enumerate(ligands):
            for ligand2 in ligands[i+1:]:
                lig_pair = LigPair(docking.ligands[ligand1],
                                   docking.ligands[ligand2],
                                   interactions, dataset.proteins[protein].lm.mcss, max_poses, True)
                ligand_stats = statistics_lig_pair(lig_pair, interactions, pnative, weighted, sd)
                protein_stats = merge_stats(protein_stats, ligand_stats)
        stats = merge_stats(stats, protein_stats)
    return stats

def gscore_statistics(data, max_poses = 100, native_thresh = 2.0, points = 1000,
                      scale_by_top = False, sd = 0.4, domain = (-16, 2)):
    '''
    data    {protein (str): [ligand (str), ...]}

    Returns DensityEstimate's representing the native and reference
    distribution of the glide scores for all proteins, ligands in data.
    '''
    stats = {}
    for protein, ligands in data.items():
        print(protein, len(ligands))
        print(ligands)
        protein_stats = {}
        dataset = Dataset(shared_paths, [protein])
        dataset.load({protein: ligands}, load_fp=False, load_mcss=False)
        st = dataset.proteins[protein].lm.st
        docking = dataset.proteins[protein].docking[st]
        for ligand in ligands:
            # Get input data
            native = np.array([pose.rmsd <= native_thresh
                for pose in docking.ligands[ligand].poses[:max_poses]])
            glide = np.array([pose.gscore
                for pose in docking.ligands[ligand].poses[:max_poses]])

            if glide.shape[0] and scale_by_top: glide -= glide.min()
            
            # Compute densities
            ligand_stats = {'native':{}, 'reference':{}}
            ligand_stats['native'][0] = DensityEstimate(points = points,
                                                        reflect = scale_by_top,
                                                        sd = sd,
                                                        domain = domain,
                                                        out_of_bounds = 0)
            ligand_stats['native'][0].fit(glide[native==1])
            ligand_stats['reference'][0] = DensityEstimate(points = points,
                                                           reflect = scale_by_top,
                                                           sd = sd,
                                                           domain = domain,
                                                           out_of_bounds = 0)
            ligand_stats['reference'][0].fit(glide)
            
            protein_stats = merge_stats(protein_stats, ligand_stats)
        stats = merge_stats(stats, protein_stats)
    return stats['native'][0], stats['reference'][0]

if __name__ == __main__:
    

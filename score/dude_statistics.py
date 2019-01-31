import numpy as np
from score.density_estimate import DensityEstimate
from score.dude_pairs import DUDEPDBLigPair
from dock.pick_helpers import load_helpers
from score.pairs import LigPair
from containers import Protein
from shared_paths import shared_paths, feature_defs

# Helper Functions
def merge_dicts_of_lists(d1, d2):
    assert (type(d1) == dict) == (type(d2) == dict)
    if type(d1) == dict:
        return {key: merge_dicts_of_lists(d1[key] if key in d1 else [],
                                          d2[key] if key in d2 else [])
                for key in set(d1.keys()).union(set(d2.keys()))}
    if type(d1) != list: d1 = [d1]
    if type(d2) != list: d2 = [d2]
    return d1 + d2

def merge_stats(stats, weight):
    for d, interactions in stats.items():
        for i, des in interactions.items():
            stats[d][i] = DensityEstimate.merge(des, weight)
    return stats

# I/O
def get_fname(protein, dude_ligand, pdb_ligand):
    '''
    Inputs:
    * protein, dude_ligand, pdb_ligand (str): names of protein and ligands
    * for protein level stats, set dude_ligand and pdb_ligand to None.

    Returns:
    * format string to paths of statistics files. The paths take the form:
        shared_paths['write_data']/protein_id/stats/statsversion/ID-{}-{}.de

        Note: if this is a protein-level statistics file, the ID is simply the protein ID. If this
        is the statistics file for a ligpair, the ID takes the form 'dude_ligand-pdb_ligand'
    '''
    assert (dude_ligand is None) == (pdb_ligand is None)
    if dude_ligand is None and pdb_ligand is None:
        ID = protein
    else:
        ID = '{}-{}'.format(dude_ligand, pdb_ligand)
    version = shared_paths['stats']['version']
    return "{}/{}/dude_stats/{}/{}-{}-{}.de".format(shared_paths['write_data'],
                                               protein, version, ID,
                                               '{}', '{}')

def read_stats(fname, interactions):
    ''' Reads all available statistics files from the correct directory and returns data in a dict.

    Inputs:
    * fname (str): output from get_fname(...)
    * interactions [interaction (str),  ...]: interactions to read from files.

    Returns:
    * stats {'dude':{'mcss': DensityEstimate object, ... },
            }
    Essentially, a dict from 'dude' to dicts mapping each feature type to a DensityEstimate
    object for that stat-type+feature-type pair
    '''
    stats = {'dude':{}}
    for stattype in stats: # stattype should only take the value 'dude', this line is sort of unnecessary
        for inttype in interactions:
            try: 
                stats[stattype][inttype] = DensityEstimate.read(fname.format(inttype, stattype))
            except:
                pass
    return stats

# Core stats computation
def get_interaction(lig_pair, interaction):
    '''
    Inputs:
    * lig_pair (LigPair): ligand pair for which to get features
    * interaction (str): interaction type to get scores for
    
    Returns:
    * X_native: features for native poses
    * X_ref:    features for all poses
    * w_ref:    pairs of glide scores for all poses
    '''
    X_native, X_ref, w_ref = [], [], []
    for (r1,r2), pp in lig_pair.pose_pairs.items():
        if max(r1, r2) > shared_paths['stats']['max_poses']: continue
        pp_x = lig_pair.get_feature(interaction, r1, r2)
        if pp_x is not None:
            if pp.correct():
                X_native += [pp_x]
            X_ref += [pp_x]

    X_native, X_ref = np.array(X_native), np.array(X_ref)
    return X_native, X_ref,  w_ref

def statistics_lig_pair(protein, dude_ligand, pdb_ligand, interactions):
    '''
    Inputs:
    * protein (str):  protein name
    * dude_ligand,2 (str): names of ligands for which to compute statistics
    * interactions [interaction (str), ...]: interactions for which to compute
        statistics
    '''
    # If all statistics for this ligpair have already been computed, just read them in and return.
    fname = get_fname(protein, dude_ligand, pdb_ligand)
    stats = read_stats(fname, interactions)
    if all(interaction in stats['dude'] for interaction in interactions):
        return stats
    
    print('Computing statistics for:', protein, dude_ligand, pdb_ligand)
    # Load relevant data.
    prot = Protein(protein)
    prot.load_docking([dude_ligand, pdb_ligand],
                      load_fp = True, load_mcss = 'mcss' in interactions)
    lm = prot.lm
    docking = prot.docking[lm.st]
    lig_pair = DUDEPDBLigPair(docking.ligands[dude_ligand],
                       docking.ligands[pdb_ligand],
                       interactions, lm.mcss if 'mcss' in interactions else None,
                       shared_paths['stats']['max_poses'])

    # Compute all remaining statistics and write to files.
    for interaction in interactions:
        if (interaction in stats['dude']):
            continue

        X_native, X_ref, w_ref = get_interaction(lig_pair, interaction)
        w_ref = 1

        for d, X in [('native', X_native), ('reference', X_ref)]:
            domain = (0, 15) if interaction == 'mcss' else (0, 1)

            stats[d][interaction] = DensityEstimate(domain = domain,
                              sd = shared_paths['stats']['stats_sd']*(domain[1]-domain[0]),
                              reflect = True)
            stats[d][interaction].fit(X)
            stats[d][interaction].write(fname.format(interaction, d))
    return stats

def statistics_protein(protein, interactions):
    # Get a format string used to produce filenames for all statistics output files
    fname = get_fname(protein, None, None)

    # Try to read any existing stats files. If they all already exist, return early!
    stats = read_stats(fname, interactions)
    if all(inttype in stats['dude'] for inttype in interactions):
        return stats
    
    print('Computing statistics for:', protein)

    # Get the first # n_ligs PDB ligands
    pdb_ligands = Protein(protein).lm.get_xdocked_ligands(shared_paths['stats']['n_ligs'])

    # Get # n_helpers DUD-E ligands
    dude_ligands = load_helpers()

    print("Computing statistics for the following PDB ligands...")
    for ligand in pdb_ligands:
        print("* {}".format(ligand))
    print("...and the following DUD-E ligands")
    for ligand in dude_ligands:
        print("* {}".format(ligand))

    # Iterate over all pairs of PDB ligands. After computing statistics for a single ligand pair,
    # merge with all previously computed stats
    for dude_ligand in dude_ligands:
        for pdb_ligand in pdb_ligands:
            ligpair_stats = statistics_lig_pair(protein, dude_ligand, pdb_ligand, interactions)
            stats = merge_dicts_of_lists(stats, ligpair_stats)

    stats = merge_stats(stats, shared_paths['stats']['poses_equal'])
    
    for d, interactions in stats.items():
        for i, de in interactions.items():
            de.write(fname.format(i, d))
    return stats

def statistics(data, interactions):
    '''
    data    {protein (str): [ligand (str), ...]}
    interactions {name (str): [code (int), ...]}

    Returns DensityEstimate's representing the native and reference
    distribution of the statistics for all proteins, ligands in data.
    '''
    stats = {'dude':{}}
    for protein in data:
        os.system("mkdir -p {}/{}/dude_stats/".format(shared_paths['write_data'], protein))
        protein_stats = statistics_protein(protein, interactions)
        stats = merge_dicts_of_lists(stats, protein_stats)
    return merge_stats(stats, shared_paths['stats']['ligands_equal'])

def main(args):
    assert len(args) == 2
    statistics([args[1]], feature_defs.keys())

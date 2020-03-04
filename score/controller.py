"""
stats_version/
              pdb/
                  standard/num-alpha-feat1_feat2_feat3
                  crystal/...
                  only_crystal/...
              best_mcss/
                        standard/
                        crystal/
              best_affinity/
                            standard/num-alpha-feat1_feat2_feat3
                            crystal/
              summary/
                  pdb.txt
                  best_mcss.txt
                  best_affinity.txt
"""

import os
from settings import proteins, feature_defs, paths, stats
from containers import Protein
from score.statistics import Statistics
from glob import glob
from utils import grouper

group_size = 1
num_ligs = [20] #[1, 3, 5, 10, 15, 20]
alpha_factors = [0.0]
features = [['mcss', 'contact', 'hbond', 'sb'],
            # ['mcss'],
            # ['mcss', 'contact', 'hbond_donor', 'hbond_acceptor', 'sb'],
            # ['mcss', 'hbond', 'sb'],
            # ['mcss', 'contact', 'hbond'],
            # ['mcss', 'contact', 'sb'],
            # ['contact', 'hbond', 'sb']
            ]

cmd = '$SCHRODINGER/run {0:}/main.py score {1:} {1:} {1:} {1:}'.format(
                                   paths['CODE'], '{}')

def score(stats_root, struct, protein, ligands, use_crystal_pose,
          alpha_factor, features, stats, crystal=False, chembl = None):
    """
    stats_root: Absolute path to statistics files.
    struct: Docking structure, likely lm.st
    protein: Name of protein
    ligands: List of PDB ligands to score
    use_crystal_pose: If true use struct's pose in scoring
    alpha_factor: Set alpha to alpha_factor * (# ligands - 1)
    features: list of features to use. Should be in feature_def.keys()
    crystal: Optimize using only the docking structure
    chembl: Optimize using chembl ligands. Set to None if don't want to use
            chembl ligands otherwise set to (pick_helpers fname, # CHEMBL).
    """
    if crystal:
        assert use_crystal_pose
        assert chembl is None
    if chembl is not None:
        assert not crystal
        assert type(chembl[0]) == str, chembl[0]
        assert type(chembl[1]) == int, chembl[1]
        assert chembl[1] > 0
    
    # Setup directory
    if chembl is not None:
        subdir = '{}-{}-{}'.format(chembl[1], alpha_factor, '_'.join(features))
    else:
        subdir = '{}-{}'.format(alpha_factor, '_'.join(features))
    if os.path.exists(subdir):
        if not glob('{}/*.sc'.format(subdir)):
            print('Missing output in {}/{}.'.format(os.getcwd(), subdir))
        return
    os.mkdir(subdir)
    os.chdir(subdir)

    settings = {}
    settings['num_poses'] = 100
    settings['stats'] = stats

    # Run specific settings.
    if chembl is not None:
        settings['num_pred_chembl'] = chembl[1]
        settings['chembl_file'] = chembl[0]
    
    settings['chembl'] = chembl is not None
    settings['alpha'] = alpha_factor
    settings['use_crystal_pose'] = use_crystal_pose
    settings['k_list'] = features
    settings['all_wrong'] = False
    settings['physics_score'] = 'gscore'
    settings['randomize'] = True
    
    with open('settings.py', 'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

    if crystal or chembl is not None:
        for i, group in enumerate(grouper(group_size, ligands)):
            with open('run{}.sh'.format(i),'w') as f:
                f.write('#!/bin/bash\n')
                for ligand in group:
                    f.write(cmd.format(stats_root, struct, protein, ligand)+'\n')
            os.system('sbatch -t 09:00:00 -p rondror run{}.sh'.format(i))
    else:
        with open('run.sh','w') as f:
            f.write('#!/bin/bash\n')
            f.write(cmd.format(stats_root, struct, protein, ' '.join(ligands))+'\n')
    
        os.system('sbatch -t 01:00:00 -p owners run.sh')
    os.chdir('..')

def compute_stats(protein, stats_root, stats, paths):
    stats_prots = [p for p in proteins if p != protein.lm.protein]
    path = paths['STATS'] + '/' + stats['version'] + '/{}-{}-{}.de'

    # Statistics computation.
    stats = Statistics(stats_prots, feature_defs, stats, paths, path)
    for dist, interactions in stats.stats.items():
        for interaction, de in interactions.items():
            fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
            if not os.path.exists(fname):
                with open(fname, 'w') as fp:
                    fp.write(str(de)+'\n')

def run_pdb(stats):
    for i, d in enumerate(proteins):
        print(d, i)
        os.chdir(paths['ROOT'].format(protein=d))
        os.system('mkdir -p {}'.format('scores'))
        
        root = '{}/{}/{}'.format(paths['ROOT'].format(protein=d), 'scores',
                                 stats['version'])
        stats_root = '{}/{}'.format(root, 'stats')
        scores_root = '{}/{}'.format(root, 'pdb')

        os.system('mkdir -p {}'.format(root))
        os.system('mkdir -p {}'.format(stats_root))
        os.system('mkdir -p {}'.format(scores_root))
        os.system('mkdir -p {}/{}'.format(scores_root, 'standard'))
        os.system('mkdir -p {}/{}'.format(scores_root, 'crystal'))
        os.system('mkdir -p {}/{}'.format(scores_root, 'only_crystal'))

        protein = Protein(d, stats, paths)
        ligands = protein.lm.get_xdocked_ligands(20)
        compute_stats(protein, stats_root, stats, paths)
        
        # Standard.
        os.chdir(scores_root)
        os.chdir('standard')
        for alpha_factor in alpha_factors:
            for feature in features:
                score(stats_root, protein.lm.st, d, ligands,
                      False, alpha_factor, feature, stats)
        os.chdir('..')

        # Crystal.
        os.chdir('crystal')
        for alpha_factor in alpha_factors:
            for feature in features:
                score(stats_root, protein.lm.st, d, ligands,
                      True, alpha_factor, feature, stats)
        os.chdir('..')

        # Crystal only.
        os.chdir('only_crystal')
        for alpha_factor in alpha_factors:
            for feature in features:
                score(stats_root, protein.lm.st, protein.lm.protein, ligands,
                      True, alpha_factor, feature, stats, crystal = True)
        os.chdir('..')

def run_chembl(stats, helpers):
    for i, d in enumerate(proteins):
        print(d, i)
        os.chdir('{}/{}'.format(paths['DATA'], d))
        os.system('mkdir -p {}'.format('scores'))
        root = '{}/{}/{}/{}'.format(paths['DATA'], d, 'scores',
                                    stats['version'])
        stats_root = '{}/{}'.format(root, 'stats')
        scores_root = '{}/{}'.format(root, helpers)

        os.system('mkdir -p {}'.format(root))
        #os.system('mkdir -p {}'.format(stats_root))
        os.system('mkdir -p {}'.format(scores_root))
        os.system('mkdir -p {}/{}'.format(scores_root, 'standard'))
        #os.system('mkdir -p {}/{}'.format(scores_root, 'crystal'))

        protein = Protein(d, stats, paths)
        ligands = protein.lm.get_xdocked_ligands(20)
        # compute_stats(protein, stats_root, stats, paths)
        # Standard.
        os.chdir(scores_root)
        os.chdir('standard')
        for num_lig in num_ligs:
            for alpha_factor in alpha_factors:
                for feature in features:
                    score(stats_root, protein.lm.st, d, ligands,
                          False, alpha_factor, feature, stats, chembl = (helpers, num_lig))
        os.chdir('..')

        # Crystal.
        #os.chdir('crystal')
        #for num_lig in num_ligs:
        #    for alpha_factor in alpha_factors:
        #        for feature in features:
        #            score(stats_root, protein.lm.st, protein.lm.protein, ligands,
        #                  True, alpha_factor, feature, stats, chembl = (helpers, num_lig))
        #os.chdir('..')

def check(stats_version, mode):

    template = '{}/*/scores/{}/{}/*/*/'.format(paths['DATA'],
                                               stats_version,
                                               mode)
    print(template)
    errors = []
    total = 0
    for fname in glob(template+'slurm*'):
        with open(fname) as fp:
            text = fp.read()
            if (('Error' in text or 'Exception' in text or 'error' in text)
                and 'PREEMPTION' not in text):
                if len(errors) < 10:
                    # Don't write too many error messages
                    print(fname)
                    print('\n'.join(text.split('\n')[-10:]))
                errors += ['/'.join(fname.split('/')[:-1])]
        total += 1
    
    print('{} slurm files'.format(total))
    print('{} scoring scripts'.format(len(glob(template+'run*.sh'))))
    print(' '.join(errors))

def merge_protein(stats_version, mode, protein):
    template = '{}/{}/scores/{}/{}/*/*/*.sc'.format(paths['DATA'],
                                                    protein,
                                                    stats_version,
                                                    mode)
    print(template)
    out_fname = '{}/{}/scores/{}/summary/{}.tsv'.format(paths['DATA'],
                                                        protein,
                                                        stats_version,
                                                        mode)

    with open(out_fname, 'w') as out:
        out.write('\t'.join(['mode', 'protein', 'ligand', 'n_ligs', 'alpha', 'features',
                             'combind_rank', 'combind_rmsd',
                             'glide_rank',   'glide_rmsd',
                             'best_rank',    'best_rmsd']) + '\n')
        for fname in glob(template):
            _protein, _, version, pdb, params, settings, name = fname.split('/')[-7:]
            settings = settings.split('-')
            if len(settings) == 2:
                settings = ['0'] + settings
            assert _protein == protein
            assert version == stats_version

            with open(fname) as fp:
                fp.readline()
                for line in fp:
                    tok = line.strip().split(',')
                    if len(tok) == 3: continue
                    (lig,
                     combind_rank, combind_rmsd,
                     glide_rank,   glide_rmsd,
                     best_rank,    best_rmsd) = tok

                    if combind_rmsd == '0': continue
                    if 'CHEMBL' in lig: continue
                    out.write('\t'.join([
                                         params,
                                         protein,
                                         lig,
                                         settings[0],
                                         settings[1],
                                         settings[2],
                                         combind_rank, combind_rmsd,
                                         glide_rank, glide_rmsd,
                                         best_rank, best_rmsd
                                        ]) + '\n')

def merge(stats_version, mode):
    for protein in proteins:
        directory = '{}/{}/scores/{}/summary'.format(paths['DATA'],
                                                     protein,
                                                     stats_version,
                                                     mode)
        print(protein)
        os.system('mkdir -p {}'.format(directory))
        merge_protein(stats_version, mode, protein)

def main(args):
    print('There are {} total PDB jobs.'.format(  len(alpha_factors)
                                                * len(features)
                                                * len(proteins)
                                                * 3))
    print('There are {} total CHEMBL jobs.'.format(  len(alpha_factors)
                                                   * len(num_ligs)
                                                   * len(features)
                                                   * len(proteins)
                                                   * 2))
    mode, stats_version = args[:2]
    name = args[2] if len(args) == 3 else None
    if   mode == 'run_pdb':
        run_pdb(stats[stats_version])
    elif mode == 'run_chembl':
        run_chembl(stats[stats_version], name)
    elif mode == 'check':
        check(stats_version, name)
    elif mode == 'merge':
        merge(stats_version, name)
    else:
        assert False

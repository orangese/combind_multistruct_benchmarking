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
from shared_paths import shared_paths, proteins, feature_defs
from containers import Protein
from score.statistics import statistics
from glob import glob
from grouper import grouper

group_size = 4
num_ligs = [1, 3, 5, 10, 15]
alpha_factors = [1.0]
features = [['mcss', 'contact', 'hbond', 'sb'],
            #['mcss'],
            ['mcss', 'hbond', 'sb'],
            #['mcss', 'contact', 'hbond'],
            #['contact', 'hbond', 'sb']
            ]

cmd = '$SCHRODINGER/run {0:}/main.py score {1:} {1:} {1:} {1:}'.format(
                                   shared_paths['code'], '{}')
settings = {
    'num_poses' : 100,
    'shared_paths': shared_paths
}

def score(stats_root, struct, protein, ligands, use_crystal_pose,
          alpha_factor, features, crystal = False, chembl = None):
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

    # Run specific settings.
    if chembl is not None:
        settings['num_pred_chembl'] = chembl[1]
        settings['chembl_file'] = chembl[0] + '.txt'
    
    settings['chembl'] = chembl is not None
    settings['alpha'] = alpha_factor
    settings['use_crystal_pose'] = use_crystal_pose
    settings['k_list'] = features
    
    
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
            os.system('sbatch -t 09:00:00 -p owners run{}.sh'.format(i))
    else:
        with open('run.sh','w') as f:
            f.write('#!/bin/bash\n')
            f.write(cmd.format(stats_root, struct, protein, ' '.join(ligands))+'\n')
    
        os.system('sbatch -t 05:00:00 -p owners run.sh')
    os.chdir('..')

def compute_stats(protein, stats_root):
    stats_prots = [p for p in proteins if p != protein.lm.protein and p != 'DHFR']

    # Statistics computation.
    stats = statistics(stats_prots, feature_defs.keys())
    for dist, interactions in stats.items():
        for interaction, de in interactions.items():
            fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
            if not os.path.exists(fname):
                with open(fname, 'w') as fp:
                    fp.write(str(de)+'\n')

def run_pdb():
    for i, d in enumerate(proteins):
        print(d, i)
        os.chdir('{}/{}'.format(shared_paths['data'], d))
        os.system('mkdir -p {}'.format('scores'))
        root = '{}/{}/{}/{}'.format(shared_paths['data'], d, 'scores',
                                    shared_paths['stats']['version'])
        stats_root = '{}/{}'.format(root, 'stats')
        scores_root = '{}/{}'.format(root, 'pdb')

        os.system('mkdir -p {}'.format(root))
        os.system('mkdir -p {}'.format(stats_root))
        os.system('mkdir -p {}'.format(scores_root))
        os.system('mkdir -p {}/{}'.format(scores_root, 'standard'))
        os.system('mkdir -p {}/{}'.format(scores_root, 'crystal'))
        os.system('mkdir -p {}/{}'.format(scores_root, 'only_crystal'))

        protein = Protein(d)
        ligands = protein.lm.get_xdocked_ligands(20)
        compute_stats(protein, stats_root)

        # Standard.
        os.chdir(scores_root)
        os.chdir('standard')
        for alpha_factor in alpha_factors:
            for feature in features:
                score(stats_root, protein.lm.st, protein.lm.protein, ligands,
                      False, alpha_factor, feature)
        os.chdir('..')

        # Crystal.
        os.chdir('crystal')
        for alpha_factor in alpha_factors:
            for feature in features:
                score(stats_root, protein.lm.st, protein.lm.protein, ligands,
                      True, alpha_factor, feature)
        os.chdir('..')

        # Crystal only.
        os.chdir('only_crystal')
        for alpha_factor in alpha_factors:
            for feature in features:
                score(stats_root, protein.lm.st, protein.lm.protein, ligands,
                      True, alpha_factor, feature, crystal = True)
        os.chdir('..')

def run_chembl(helpers):
    for i, d in enumerate(proteins):
        print(d, i)
        os.chdir('{}/{}'.format(shared_paths['data'], d))
        os.system('mkdir -p {}'.format('scores'))
        root = '{}/{}/{}/{}'.format(shared_paths['data'], d, 'scores',
                                 shared_paths['stats']['version'])
        stats_root = '{}/{}'.format(root, 'stats')
        scores_root = '{}/{}'.format(root, helpers)

        os.system('mkdir -p {}'.format(root))
        os.system('mkdir -p {}'.format(stats_root))
        os.system('mkdir -p {}'.format(scores_root))
        os.system('mkdir -p {}/{}'.format(scores_root, 'standard'))
        os.system('mkdir -p {}/{}'.format(scores_root, 'crystal'))

        protein = Protein(d)
        ligands = protein.lm.get_xdocked_ligands(20)
        compute_stats(protein, stats_root)
        # Standard.
        os.chdir(scores_root)
        os.chdir('standard')
        for num_lig in num_ligs:
            for alpha_factor in alpha_factors:
                for feature in features:
                    score(stats_root, protein.lm.st, protein.lm.protein, ligands,
                          False, alpha_factor, feature, chembl = (helpers, num_lig))
        os.chdir('..')

        # Crystal.
        os.chdir('crystal')
        for num_lig in num_ligs:
            for alpha_factor in alpha_factors:
                for feature in features:
                    score(stats_root, protein.lm.st, protein.lm.protein, ligands,
                          True, alpha_factor, feature, chembl = (helpers, num_lig))
        os.chdir('..')

def check(mode):

    template = '{}/*/scores/{}/{}/*/*/'.format(shared_paths['data'],
                                               shared_paths['stats']['version'],
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

def merge_protein(mode, protein):
    """
    Unlike the other functions here merge acts on only one protein.
    Doing all proteins at once simply takes too long.
    """
    template = '{}/{}/scores/{}/{}/*/*/*.sc'.format(shared_paths['data'],
                                                    protein,
                                                    shared_paths['stats']['version'],
                                                    mode)

    print(template)

    out_fname = '{}/{}/scores/{}/summary/{}.tsv'.format(shared_paths['data'],
                                                        protein,
                                                        shared_paths['stats']['version'],
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
            assert version == shared_paths['stats']['version']
                
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

def merge(mode, inline):
    for protein in proteins:
        directory = '{}/{}/scores/{}/summary'.format(shared_paths['data'],
                                                     protein,
                                                     shared_paths['stats']['version'],
                                                     mode)
        log = '{}/{}.log'.format(directory, mode)
        print(protein)
        if os.path.exists(log):
            print('{} exists:'.format(log))
            with open(log) as fp:
                print(fp.read())
            continue
        print('executing:')
        os.system('mkdir -p {}'.format(directory))
        if inline:
            merge_protein(mode, protein)
            with open(log, 'w') as fp:
                pass
        else:
            os.system('sbatch -p owners -t 01:00:00 --out={} '
                      '--wrap="python {}/score/score_controller.py '
                      'merge_protein {} {}"'.format(log, shared_paths['code'], mode, protein))

def archive(mode):
    template = '{}/*/scores/{}/{}'.format(shared_paths['data'],
                                          shared_paths['stats']['version'],
                                          mode)
    paths = glob(template)
    for path in paths:
        print(path)
    if input('Archive all of the above directories? (yes)') != 'yes':
        return

    for path in glob(template):
        os.system('tar -zcvf {0:}.tar.gz {0:}'.format(path))
        os.system('rm -rd {}'.format(path))

def inflate(mode):
    template = '{}/*/scores/{}/{}.tar.gz'.format(shared_paths['data'],
                                                 shared_paths['stats']['version'],
                                                 mode)
    paths = glob(template)
    for path in paths:
        print(path)
    if input('Inflate all of the above tarballs? (yes)') != 'yes':
        return

    for path in glob(template):
        if os.path.exists(path.replace('.tar.gz', '')):
            print(path)
            continue
        os.system('tar -xzf {}'.format(path))
        os.system('rm {}'.format(path))

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
    if   args[1] == 'run_pdb':
        run_pdb()
    elif args[1] == 'run_chembl':
        run_chembl(args[2])
    elif args[1] == 'check':
        check(args[2])
    elif args[1] == 'merge':
        merge(args[2], len(args) == 4 and args[3] == 'inline')
    elif args[1] == 'merge_protein':
        merge_protein(args[2], args[3])
    elif args[1] == 'archive':
        archive(args[2])
    elif args[1] == 'inflate':
        inflate(args[2])
    else:
        assert False

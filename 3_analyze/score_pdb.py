"""
stats_version/
              pdb/
                  standard/alpha-feat1_feat2_feat3
                  crystal/...
                  only_crystal/...
"""

import os
from shared_paths import shared_paths, proteins, feature_defs
from containers import Protein
from statistics import statistics
from glob import glob

alpha_factors = [0.5, 1.0, 1.5, 2.0, 2.5]
features = [['mcss', 'pipi', 'contact', 'hbond', 'sb'],
            ['mcss', 'contact', 'hbond', 'sb'],
            ['pipi', 'contact', 'hbond', 'sb'],
            ['mcss', 'pipi', 'contact', 'sb'],
            ['mcss', 'pipi', 'hbond', 'sb'],
            ['mcss', 'pipi', 'contact', 'hbond'],
            ]

cmd = '$SCHRODINGER/run {0:}/3_analyze/scores.py {1:} {1:} {1:} {1:}'.format(
                                                  shared_paths['code'], '{}')
settings = {
    'num_poses' : 100,
    'chembl': False,
    'shared_paths': shared_paths
}

def score(stats_root, struct, protein, ligands, use_crystal_pose,
          alpha_factor, features, crystal = False):
    subdir = '{}-{}'.format(alpha_factor, '_'.join(features))
    if os.path.exists(subdir):
        if not glob('{}/*.sc'.format(subdir)):
            print('Missing output in {}/{}.'.format(os.getcwd(), subdir))
        return
    os.mkdir(subdir)
    os.chdir(subdir)
    
    if crystal:
        settings['alpha'] = alpha_factor
    else:
        settings['alpha'] = float(len(ligands)-1) * alpha_factor
    settings['use_crystal_pose'] = use_crystal_pose
    settings['k_list'] = features


    
    with open('settings.py', 'w') as f:
        for varname, var in settings.items():
            if type(var) is str: var = '"{}"'.format(var)
            f.write('{}={}\n'.format(varname, var))

    with open('run.sh','w') as f:
        f.write('#!/bin/bash\n')
        if crystal:
            assert use_crystal_pose
            for ligand in ligands:
                f.write(cmd.format(stats_root, struct, protein, ligand)+'\n')
        else:
            f.write(cmd.format(stats_root, struct, protein, ' '.join(ligands))+'\n')
    os.system('sbatch -t 01:00:00 -p owners run.sh')
    os.chdir('..')

print('There are {} total jobs.'.format(  len(alpha_factors)
                                        * len(features)
                                        * len(proteins)
                                        * 3))

def run():
    for i, d in enumerate(proteins):
        print(d, i)
        os.chdir('{}/{}'.format(shared_paths['data'], d))
        os.system('mkdir -p {}'.format('scores'))
        root = '{}/{}/{}/{}'.format(shared_paths['data'], d, 'scores',
                                 shared_paths['stats']['version'])
        stats_root = '{}/{}'.format(root, 'stats')
        pdb_root = '{}/{}'.format(root, 'pdb')

        os.system('mkdir -p {}'.format(root))
        os.system('mkdir -p {}'.format(stats_root))
        os.system('mkdir -p {}'.format(pdb_root))
        os.system('mkdir -p {}/{}'.format(pdb_root, 'standard'))
        os.system('mkdir -p {}/{}'.format(pdb_root, 'crystal'))
        os.system('mkdir -p {}/{}'.format(pdb_root, 'only_crystal'))

        protein = Protein(d)
        ligands = protein.lm.get_xdocked_ligands(20)
        stats_prots = [p for p in proteins if p != protein.lm.protein]

        # Statistics computation.
        stats = statistics(stats_prots, feature_defs.keys())
        for dist, interactions in stats.items():
            for interaction, de in interactions.items():
                with open('{}/{}_{}.txt'.format(stats_root, dist, interaction), 'w') as fp:
                    fp.write(str(de)+'\n')

        # Standard.
        os.chdir(pdb_root)
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
                      True, alpha_factor, feature, True)
        os.chdir('..')

def check():

    template = '{}/*/scores/{}/*/*/*/'.format(shared_paths['data'],
                                                  shared_paths['stats']['version'])
    print(template)
    errors = []
    total = 0
    for fname in glob(template+'slurm*'):
        with open(fname) as fp:
            text = fp.read()
            if 'Error' in text or 'Exception' in text or 'error' in text:
                if len(errors) < 10:
                    # Don't write too many error messages
                    print(fname)
                    print(text)
                errors += ['/'.join(fname.split('/')[:-1])]
        total += 1
    
    print('{} slurm files'.format(total))
    print('{} scoring scripts'.format(len(glob(template+'run.sh'))))
    print(' '.join(errors))

def merge():
    template = '{}/*/scores/{}/pdb/*/*/*.sc'.format(shared_paths['data'],
                                                    shared_paths['stats']['version'])

    out_fname = '{}/../bpp_outputs/{}.tsv'.format(shared_paths['data'],
                                                  shared_paths['stats']['version'])
    total = 0
    with open(out_fname, 'w') as out:
        out.write('\t'.join(['mode', 'protein', 'ligand', 'alpha', 'features',
                             'combind_rank', 'combind_rmsd',
                             'glide_rank',   'glide_rmsd',
                             'best_rank',    'best_rmsd']) + '\n')
        for fname in glob(template):
            protein, _, version, pdb, mode, settings, name = fname.split('/')[-7:]
            settings = settings.split('-')
            assert len(settings) == 2
            assert pdb == 'pdb'
            assert version == shared_paths['stats']['version']
            with open(fname) as fp:
                fp.readline()
                for line in fp:
                    tok = line.strip().split(',')
                    if len(tok) == 3:
                        continue
                    else:
                        (lig,
                         combind_rank, combind_rmsd,
                         glide_rank,   glide_rmsd,
                         best_rank,    best_rmsd) = tok
                        if combind_rmsd == '0': continue
                        if 'CHEMBL' in lig: continue
                        out.write('\t'.join([
                                            mode,
                                            protein,
                                            lig,
                                            settings[0],
                                            settings[1],
                                            combind_rank, combind_rmsd,
                                            glide_rank, glide_rmsd,
                                            best_rank, best_rmsd
                                            ]) + '\n')
                        total += 1
    print('Processed a total of {} entries.'.format(total))

def archive():
    template = '{}/*/scores/{}'.format(shared_paths['data'],
                                       shared_paths['stats']['version'])
    paths = glob(template)
    for path in paths:
        print(path)
    if input('Archive all of the above directories? (yes)') != 'yes':
        return

    for path in glob(template):
        if os.path.exists(path + '.tar.gz'):
            print(path)
            continue
        os.system('tar -zcvf {0:}.tar.gz {0:}'.format(path))
        os.system('rm -rd {}'.format(path))

def inflate():
    template = '{}/*/scores/{}.tar.gz'.format(shared_paths['data'],
                                              shared_paths['stats']['version'])
    paths = glob(template)
    for path in paths:
        print(path)
    if input('Inflate all of the above tarballs? (yes)') != 'yes':
        return

    for path in glob(template):
        if os.path.exists(path.replace('.tar.gz', '')):
            print(path)
            continue
        os.system('tar -xvzf {}.tar.gz'.format(path))
        os.system('rm {}'.format(path))

if __name__ ==  '__main__':
    import sys
    if   sys.argv[1] == 'run':
        run()
    elif sys.argv[1] == 'check':
        check()
    elif sys.argv[1] == 'merge':
        merge()
    elif sys.argv[1] == 'archive':
        archive()
    elif sys.argv[1] == 'inflate':
        inflate()
    else:
        assert False

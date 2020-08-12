"""
stats_version/
              pdb/
                  standard/num-alpha-feat1_feat2_feat3
                  crystal/...
                  only_crystal/...
              best_affinity/
                        standard/num-alpha-feat1_feat2_feat3
                        crystal/
              summary/
                  pdb.txt
                  best_mcss.txt
                  best_affinity.txt

for i in $(ls --color=none /oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems); do python scripts/benchmark_pose_pred.py stats rd1 /oak/stanford/groups/rondror/users/jpaggi/temp/shape $i; done;
for i in $(ls --color=none /oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems); do python scripts/benchmark_pose_pred.py setup-pdb rd1  $i --features shape; done;
n02 login ~/combind]$ for i in /oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems/*/scores/rd1/pdb/standard/0.5-shape/*.sh; do if [ ! -f ${i/.sh/.sc} ]; then cd ${i%/*}; sh ${i##*/}; cd /oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems ; fi; done;
python scripts/performance.py /oak/stanford/groups/rondror/users/jpaggi/combind/*/scores/rd2/pdb/*/*/*.sc
"""

import os
import sys
sys.path.append(os.environ['COMBINDHOME'])

import numpy as np
import pandas as pd
import click
import subprocess
from glob import glob
import config
import utils

@click.group()
def main():
    pass

# python scripts/benchmark_pose_pred.py stats rd1 /oak/stanford/groups/rondror/users/jpaggi/temp/rd1 B1AR
@main.command()
@click.argument('stats_version')
@click.argument('stats_root')
@click.argument('protein')
@click.option('--data', default='/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def stats(stats_version, stats_root, protein, data, ligands):
    import score.statistics

    merged_root = '{data}/{protein}/scores/{stats_version}/stats'
    merged_root = merged_root.format(data=data, protein=protein,
                                     stats_version=stats_version)

    paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
             'DATA': data,
             'PDB': ligands}
    paths.update(config.PATHS)
    paths = utils.resolve(paths)
    params = config.STATS[stats_version]
    
    proteins = utils.get_proteins(paths, [protein])

    score.statistics.compute(params, paths, config.FEATURE_DEFS, stats_root,
                             proteins, merged_root)

# python scripts/benchmark_pose_pred.py setup-pdb rd1 B1AR
@main.command()
@click.argument('stats_version')
@click.argument('protein')
@click.option('--alpha', default=1.0)
@click.option('--features', default='mcss,hbond,sb,contact')
@click.option('--data', default='/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def setup_pdb(stats_version, protein, alpha, features, data, ligands):
    cwd = '{data}/{protein}/scores/{stats_version}/pdb/standard/{alpha}-{features}'
    cmd = ('python /home/users/jpaggi/combind/main.py '
           '--ligands {ligands} --data {data} score {protein} cross '
           '--stats-version {stats_version} '
           '--stats-root {data}/{protein}/scores/{stats_version}/stats '
           '--alpha {alpha} --features {features} '
           '--pose-fname poses.sc  --num-poses 100')

    cwd = cwd.format(data=data, protein=protein, stats_version=stats_version,
                     alpha=alpha, features=features.replace(',', '_'))

    cmd = cmd.format(data=data, ligands=ligands, protein=protein,
                     stats_version=stats_version, alpha=alpha, features=features)

    os.makedirs(cwd, exist_ok=True)
    with open(cwd + '/poses.sh', 'w') as fp:
        fp.write(cmd + '\n')

# python scripts/benchmark_pose_pred.py setup-pdb-xtal rd1 B1AR
@main.command()
@click.argument('stats_version')
@click.argument('protein')
@click.option('--alpha', default=1.0)
@click.option('--features', default='mcss,hbond,sb,contact')
@click.option('--data', default='/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind_paper_systems')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def setup_pdb_xtal(stats_version, protein, alpha, features, data, ligands):
    cwd = '{data}/{protein}/scores/{stats_version}/pdb/crystal/{alpha}-{features}'
    cmd = ('python /home/users/jpaggi/combind/main.py '
           '--ligands {ligands} --data {data} score {protein} all '
           '--stats-version {stats_version} '
           '--stats-root {data}/{protein}/scores/{stats_version}/stats '
           '--alpha {alpha} --features {features} '
           '--pose-fname poses.sc  --num-poses 100 --xtal {st}_lig')

    st = glob('{data}/{protein}/docking/grids/*'.format(data=data, protein=protein))
    print('{data}/{protein}/docking/grids/*'.format(data=data, protein=protein))
    assert len(st) == 1, st
    st = st[0].split('/')[-1]

    cwd = cwd.format(data=data, protein=protein, stats_version=stats_version,
                     alpha=alpha, features=features.replace(',', '_'))

    cmd = cmd.format(data=data, ligands=ligands, protein=protein, st=st,
                     stats_version=stats_version, alpha=alpha, features=features)

    os.makedirs(cwd, exist_ok=True)
    with open(cwd + '/poses.sh', 'w') as fp:
        fp.write(cmd + '\n')

@main.command()
@click.argument('stats_version')
@click.argument('protein')
@click.option('--alpha', default=1.0)
@click.option('--features', default='mcss,hbond,sb,contact')
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/combind')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def setup_pdb_xtal_only(stats_version, protein, alpha, features, data, ligands):
    cwd = '{data}/{protein}/scores/{stats_version}/pdb/crystal_only/{alpha}-{features}'
    cmd = ('python /home/users/jpaggi/combind/main.py '
           '--ligands {ligands} --data {data} score {protein} {st}_lig {lig}_lig '
           '--stats-version {stats_version} '
           '--stats-root {data}/{protein}/scores/{stats_version}/stats '
           '--alpha {alpha} --features {features} '
           '--pose-fname {lig}.sc  --num-poses 100 --xtal {st}_lig')

    st = glob('{data}/{protein}/docking/grids/*'.format(data=data, protein=protein))
    assert len(st) == 1, st
    st = st[0].split('/')[-1]

    cwd = cwd.format(data=data, protein=protein, stats_version=stats_version,
                     alpha=alpha, features=features.replace(',', '_'))

    ligs = pd.read_csv(ligands.format(ROOT='{}/{}'.format(data, protein)))
    for lig in ligs["ID"]:
        if lig == st: continue
        
        _cmd = cmd.format(data=data, ligands=ligands, protein=protein, st=st, lig=lig,
                          stats_version=stats_version, alpha=alpha, features=features)

        os.makedirs(cwd, exist_ok=True)
        with open('{}/{}.sh'.format(cwd, lig), 'w') as fp:
            fp.write(_cmd + '\n')

@main.command()
@click.argument('stats_version')
@click.argument('protein')
@click.option('--alpha', default=1.0)
@click.option('--features', default='mcss,hbond,sb,contact')
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/combind')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def setup_pdb_xtal_all(stats_version, protein, alpha, features, data, ligands):
    cwd = '{data}/{protein}/scores/{stats_version}/pdb/crystal_all/{alpha}-{features}'
    cmd = ('python /home/users/jpaggi/combind/main.py '
           '--ligands {ligands} --data {data} score {protein} all '
           '--stats-version {stats_version} '
           '--stats-root {data}/{protein}/scores/{stats_version}/stats '
           '--alpha {alpha} --features {features} '
           '--pose-fname {lig}.sc  --num-poses 100 {xtal}')

    st = glob('{data}/{protein}/docking/grids/*'.format(data=data, protein=protein))
    assert len(st) == 1, st
    st = st[0].split('/')[-1]

    cwd = cwd.format(data=data, protein=protein, stats_version=stats_version,
                     alpha=alpha, features=features.replace(',', '_'))

    ligs = pd.read_csv(ligands.format(ROOT='{}/{}'.format(data, protein)))
    for lig in ligs["ID"]:
        if lig == st: continue

        xtal = list(ligs["ID"])
        xtal.remove(lig)
        xtal = ''.join('--xtal {}_lig '.format(_xtal) for _xtal in xtal)
        
        _cmd = cmd.format(data=data, ligands=ligands, protein=protein, st=st, lig=lig,
                          stats_version=stats_version, alpha=alpha, features=features,
                          xtal=xtal)

        os.makedirs(cwd, exist_ok=True)
        with open('{}/{}.sh'.format(cwd, lig), 'w') as fp:
            fp.write(_cmd + '\n')

# python scripts/benchmark_pose_pred.py setup-chembl affinity_diverse rd1 B1AR
@main.command()
@click.argument('mode')
@click.argument('stats_version')
@click.argument('protein')
@click.option('--n-helpers', default=20)
@click.option('--alpha', default=1.0)
@click.option('--features', default='mcss,hbond,sb,contact')
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/combind')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def setup_chembl(mode, stats_version, protein, n_helpers, alpha, features, data, ligands):
    cwd = '{data}/{protein}/scores/{stats_version}/{mode}/standard/{n_helpers}-{alpha}-{features}'
    cmd = ('python /home/users/jpaggi/combind/main.py '
           '--ligands {query}.csv --data {data} score {protein} all '
           '--stats-version {stats_version} '
           '--stats-root {data}/{protein}/scores/{stats_version}/stats '
           '--alpha {alpha} --features {features} '
           '--pose-fname {query}.sc --num-poses 100')

    helpers = '{data}/{protein}/chembl/*-{mode}.csv'
    helpers = helpers.format(data=data, protein=protein, mode=mode)

    for _helpers in glob(helpers):
        query = _helpers.split('/')[-1].split('-')[0]

        _cwd = cwd.format(mode=mode, data=data, protein=protein, n_helpers=n_helpers,
                          stats_version=stats_version, alpha=alpha,
                          features=features.replace(',', '_'))

        _cmd = cmd.format(data=data, helpers=_helpers, protein=protein,
                  stats_version=stats_version, alpha=alpha, query=query,
                  features=features)
        
        os.makedirs(_cwd, exist_ok=True)

        # Write version of helpers with specified number of helpers
        np.random.seed(hash(protein)%(2**32 - 1))
        helpers_df = pd.read_csv(_helpers)
        query_mask = helpers_df['ID'] == query
        if n_helpers < len(helpers_df)-1:
            helpers_df = pd.concat([helpers_df[query_mask],
                                    helpers_df[~query_mask].sample(n_helpers)])
        helpers_df.to_csv('{}/{}.csv'.format(_cwd, query), index=False)
        
        with open('{}/{}.sh'.format(_cwd, query), 'w') as fp:
            fp.write(_cmd + '\n')

@main.command()
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/combind')
@click.argument('mode')
@click.argument('stats_version')
@click.option('--scoring', default='standard')
@click.option('--inline', is_flag=True)
def run(mode, stats_version, data, scoring, inline):
    shs = '{data}/*/scores/{stats_version}/{mode}/{scoring}/*/*.sh'
    shs = shs.format(data=data, stats_version=stats_version, mode=mode, scoring=scoring)
    for sh in glob(shs):
        cwd = '/'.join(sh.split('/')[:-1])
        name = sh.split('/')[-1].split('.')[0]

        cmd = ('sbatch -p owners -t 01:00:00 -D {cwd} -J {name} '
               '-o {name}.log --wrap="sh {name}.sh"').format(name=name, cwd=cwd)
        
        if not os.path.exists('{}/{}.sc'.format(cwd, name)):
            print(cmd)
            if inline:
                subprocess.call('sh {name}.sh > {name}.log'.format(name=name), cwd=cwd, shell=True)
            else:
                os.system(cmd)

@main.command()
@click.argument('stats_version')
@click.argument('mode')
@click.option('--data', default='/oak/stanford/groups/rondror/users/jpaggi/combind')
@click.option('--ligands', default='{ROOT}/structures/pdb.csv')
def merge(stats_version, mode, data, ligands):
    paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
             'DATA': data,
             'PDB': ligands}
    paths.update(config.PATHS)
    paths = utils.resolve(paths)
    params = config.STATS[stats_version]
    
    proteins = utils.get_proteins(paths, [])
    for protein in proteins:
        directory = '{}/{}/scores/{}/summary'.format(data,
                                                     protein,
                                                     stats_version,
                                                     mode)
        print(protein)
        os.system('mkdir -p {}'.format(directory))
        merge_protein(stats_version, mode, protein, data)

def merge_protein(stats_version, mode, protein, data):
    template = '{}/{}/scores/{}/{}/*/*/*.sc'.format(data,
                                                    protein,
                                                    stats_version,
                                                    mode)
    print(template)
    out_fname = '{}/{}/scores/{}/summary/{}.tsv'.format(data,
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

            ligands = {}
            with open(fname) as fp:
                fp.readline()
                for line in fp:
                    tok = line.strip().split(',')
                    if len(tok) == 3: continue
                    if tok[2] == 'None': continue
                    ligands[tok[0]] = tok[1:]

            query = name.replace('.sc', '_lig')
            if query in ligands:
                ligands = {query: ligands[query]}

            for lig, tok in ligands.items():
                (combind_rank, combind_rmsd,
                 glide_rank,   glide_rmsd,
                 best_rank,    best_rmsd) = tok
                out.write('\t'.join([params, protein, lig, settings[0],
                                     settings[1], settings[2],
                                     combind_rank, combind_rmsd,
                                     glide_rank, glide_rmsd,
                                     best_rank, best_rmsd]) + '\n')

main()

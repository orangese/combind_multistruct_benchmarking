"""
for i in *; do python ~/combind/dock/ifd.py setup $i; done;
for i in */docking/ifd/*; do cd $i; if [ ! -f *.log ]; then echo $i; sbatch -p owners -t 12:00:00 -n 6 --wrap="sh ${i##*/}.sh" -J ${i##*/}; fi; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;
for i in */docking/ifd/*; do cd $i; echo $i;  if [ -f ${i##*/}-out.maegz ]; then python ~/combind/dock/ifd.py extract *_out.mae ${i##*/}-out.maegz ${i##*/}_pv.maegz; fi; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;
for i in */docking/ifd/*; do cd $i; echo $i;  if [ -f ${i##*/}-out.maegz ]; then python ~/combind/dock/ifd.py rmsd ${i##*/}_pv.maegz; fi; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;

# start to end of glide1
# check that all glide1 subjobs complete
# prime (in serial) through end of glide2
# check that all glide2 subjobs complete
# score.

"""

import os
from datetime import datetime
import click
import subprocess
from glob import glob
import pandas as pd


template="""INPUT_FILE  {grid}_out.mae

STAGE VDW_SCALING
  BINDING_SITE ligand L:{resid}

STAGE PREDICT_FLEXIBILITY
  BINDING_SITE ligand L:{resid}

STAGE INITIAL_DOCKING
  BINDING_SITE ligand L:{resid}
  INNERBOX 15.0
  OUTERBOX 30.0
  LIGAND_FILE {ligand}.mae
  LIGANDS_TO_DOCK all
  VARIANTS_TO_RUN A,B,C,D,E,F,G
  DOCKING_RINGCONFCUT 2.5
  DOCKING_AMIDE_MODE penal

STAGE COMPILE_RESIDUE_LIST
  DISTANCE_CUTOFF 5.0

STAGE PRIME_REFINEMENT
  NUMBER_OF_PASSES  1
  USE_MEMBRANE no
  OPLS_VERSION OPLS3e

STAGE GLIDE_DOCKING2
  BINDING_SITE ligand Z:999
  INNERBOX 5.0
  OUTERBOX auto
  LIGAND_FILE {ligand}.mae
  LIGANDS_TO_DOCK existing
  DOCKING_PRECISION SP
  DOCKING_CANONICALIZE False
  DOCKING_RINGCONFCUT 2.5
  DOCKING_AMIDE_MODE penal

STAGE SCORING
  SCORE_NAME  r_psp_IFDScore
  TERM 1.000,r_psp_Prime_Energy,1
  TERM 9.057,r_i_glide_gscore,0
  TERM 1.428,r_i_glide_ecoul,0
  REPORT_FILE report.csv
"""

def ifd(ligand, grid, root):
    from schrodinger.structure import StructureReader
    name = '{}-to-{}'.format(ligand, grid)
    cwd = '{}/docking/ifd2/{}'.format(root, name)

    cmd = ('{}/ifd -NGLIDECPU 1 -NPRIMECPU 6 {}.inp -HOST localhost'
           ' -SUBHOST localhost -WAIT -STRICT -RESTART'
           ).format(os.environ['SCHRODINGER'], name)

    if os.path.exists(cwd):
        print('{} exists. Not overwriting.'.format(cwd))
        return

    with StructureReader('{}/structures/ligands/{}_lig.mae'.format(root, grid)) as st:
        st = list(st)
        assert len(st) == 1
        st = st[0]

        residue = list(st.residue)
        assert len(residue) == 1
        residue = residue[0]

        assert residue.chain == 'L'
        resid = residue.resnum

    os.mkdir(cwd)

    with StructureReader('{}/ligands/{}/{}.mae'.format(root, ligand, ligand)) as st:
        st = list(st)[0]
        st.write('{}/{}.mae'.format(cwd, ligand))

    os.system('cp {}/structures/aligned/{}/rot-{}_query.mae {}/{}_out.mae'.format(root, grid, grid, cwd, grid))

    with open('{}/{}.inp'.format(cwd, name), 'w') as fp:
        fp.write(template.format(grid=grid, ligand=ligand, resid=resid))

    with open('{}/{}.sh'.format(cwd, name), 'w') as fp:
        fp.write(cmd + '\n')

def get_all_jobs(pattern):
    jobs = []
    for job in glob(pattern):
        job = job.split('/')[-1]
        jobs += [job]
    return set(jobs)

def get_running_jobs(all_jobs, running=True):
    cmd = ['squeue', '-u', os.environ['USER'], '-o', '%.20j']
    if running:
        cmd += ['-t', 'RUNNING']
    slurm = subprocess.run(cmd, capture_output=True, encoding='utf-8')

    jobs = []
    for job in slurm.stdout.split('\n'):
        job = job.strip()
        if job in all_jobs:
            jobs += [job]
    return set(jobs)

def get_completed_jobs(pattern):
    pattern += '/*-out.maegz'
    jobs = []
    for job in glob(pattern):
        jobs += [job.split('/')[-1].split('-out.maegz')[0]]
    return set(jobs)

@click.group()
def main():
    pass

def wildcard(path):
    resolved = glob(path)
    assert len(resolved) < 2, (path, resolved)

    if not resolved:
        return None
    return resolved[0]

@main.command()
@click.argument('workdir')
def check_initial(workdir):
    failed = []
    for i, a in enumerate(['A','B','C','D','E','F','G']):
        pv  = wildcard('{}/*scale_lig1_G_batchglide_0000{}_pv.maegz'.format(workdir, i))
        log = wildcard('{}/*scale_lig1_{}.log'.format(workdir, a))

        if pv:
            continue

        with open(log) as fp:
            txt = fp.read()

        phrases = ['** NO ACCEPTABLE LIGAND POSES WERE FOUND **',
                   'NO VALID POSES AFTER MINIMIZATION: SKIPPING.']

        if any(phrase in txt for phrase in phrases):
            continue

        failed += [a]

    print(workdir, failed)

@main.command()
@click.argument('root')
def setup(root):
    from schrodinger.structure import StructureReader
    df = pd.read_csv('{}/structures/pdb.csv'.format(root))
    grid = [fname for fname in os.listdir('{}/docking/grids'.format(root)) if fname[0] != '.']
    assert len(grid) ==  1
    grid = grid[0]
    print(grid)
    for _, row in df.iterrows():
        print(row['ID'])
        ifd(row['ID']+'_lig', grid, root)

@main.command()
@click.argument('pattern')
def run(pattern):
    n_jobs = 100
    print('Launching {} at {}'.format(pattern, datetime.now()))

    all_jobs = get_all_jobs(pattern)
    queued_jobs = get_running_jobs(all_jobs, running=False)
    running_jobs = get_running_jobs(all_jobs)
    completed_jobs = get_completed_jobs(pattern)
    remaining_jobs = all_jobs.difference(queued_jobs.union(completed_jobs))

    print('{} total jobs.'.format(len(all_jobs)))
    print('{} completed jobs.'.format(len(completed_jobs)))
    print('{} queued jobs.'.format(len(queued_jobs)))
    print('{} running jobs.'.format(len(running_jobs)))
    print('{} remaining jobs.'.format(len(remaining_jobs)))

    to_submit = min(len(remaining_jobs), max(0, n_jobs-len(queued_jobs)))
    print('Submitting {} more jobs.'.format(to_submit))

    for job in sorted(remaining_jobs)[:to_submit]:
        cmd = 'sbatch -p owners -t 12:00:00 -n 6 --wrap="sh {0}.sh" -J {0}'.format(job)
        cwd = glob(pattern[:-1]+job)[0]

        print(cwd)
        print(cmd)

        subprocess.run(cmd, cwd=cwd, shell=True)

@main.command()
@click.argument('receptor')
@click.argument('ifd-out')
@click.argument('pv')
def extract(receptor, ifd_out, pv):
    from schrodinger.structure import StructureReader
    with StructureReader(receptor) as st:
        st = list(st)
        assert len(st) == 1
        st = st[0]

        idx = None
        for chain in st.chain:
            if chain.name == 'L':
                idx = chain.getAtomIndices()
                break
        if idx:
            st.deleteAtoms(idx)

        st.write(pv)

    with StructureReader(ifd_out) as sts:
        for st in sts:
            for chain in st.chain:
                if chain.name == 'Z':
                    idx = chain.getAtomIndices()
                    st = st.extract(idx)
                    st.append(pv)
                    break
            else:
                assert False

@main.command()
@click.argument('pv')
def rmsd(pv):
    ligand = pv.split('-to-')[0]
    os.system('run rmsd.py -use_neutral_scaffold -pv second -c rmsd.csv '
              '../../../structures/ligands/{}.mae {}'.format(ligand, pv))

main()

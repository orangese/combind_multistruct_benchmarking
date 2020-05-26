"""
for i in *; do python ~/combind/dock/ifd.py setup $i; done;
for i in */docking/ifd/*; do cd $i; if [ ! -f *.log ]; then echo $i; sbatch -p owners -t 12:00:00 -n 6 --wrap="sh ${i##*/}.sh" -J ${i##*/}; fi; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;
for i in */docking/ifd/*; do cd $i; echo $i;  if [ -f ${i##*/}-out.maegz ]; then python ~/combind/dock/ifd.py extract *_out.mae ${i##*/}-out.maegz ${i##*/}_pv.maegz; fi; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;
for i in */docking/ifd/*; do cd $i; echo $i;  if [ -f ${i##*/}-out.maegz ]; then python ~/combind/dock/ifd.py rmsd ${i##*/}_pv.maegz; fi; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;

# Schrodinger interpreters
for i in *; do python ~/combind/dock/ifd2.py setup-stage1 $i; done;
for i in *; do python ~/combind/dock/ifd2.py setup-stage2 /oak/stanford/groups/rondror/users/jpaggi/combind/$i/docking/ifd2; done;
for i in *; do python ~/combind/dock/ifd2.py setup-stage3 /oak/stanford/groups/rondror/users/jpaggi/combind/$i/docking/ifd2; done;


# python 3.8 interpreter (conda activate mol)
python ~/combind/dock/ifd2.py run '*/docking/ifd2/*/*-stage1.inp' '*/docking/ifd2/*/*-stage1.log'
python ~/combind/dock/ifd2.py run '*/docking/ifd2/*/*-stage2.inp' '*/docking/ifd2/*/*-stage2.log' --time 12:00:00


for i in */docking/ifd2/*; do cd $i; echo $i;  python ~/combind/dock/ifd2.py extract *_out.mae ${i##*/}-stage1-out.maegz ${i##*/}_pv.maegz; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;
for i in */docking/ifd2/*; do cd $i; echo $i;  python ~/combind/dock/ifd2.py rmsd ${i##*/}_pv.maegz; cd /oak/stanford/groups/rondror/users/jpaggi/combind; done;

"""

import os
from datetime import datetime
import click
import subprocess
from glob import glob
import pandas as pd

def wildcard(path):
    resolved = glob(path)
    assert len(resolved) < 2, (path, resolved)

    if not resolved:
        return None
    return resolved[0]

def get_all_jobs(pattern):
    jobs = {}
    for job in glob(pattern):
        name = job.split('/')[-1].split('.')[0]
        cwd = '/'.join(job.split('/')[:-1])
        jobs[name] = cwd
    return jobs

def get_running_jobs(all_jobs, running=True):
    cmd = ['squeue', '-u', os.environ['USER'], '-o', '%.50j']
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
    pattern += '/*-stage1.log'
    return get_all_jobs(pattern)

@click.group()
def main():
    pass

@main.command()
@click.argument('input_pattern')
@click.argument('completed_pattern')
@click.option('--time', default='2:00:00')
@click.option('--queue', default='owners')
@click.option('--n-jobs', default=200)
def run(input_pattern, completed_pattern, time, queue, n_jobs):
    
    print('Launching {} at {}'.format(input_pattern, datetime.now()))

    all_jobs = get_all_jobs(input_pattern)
    
    queued_jobs = get_running_jobs(all_jobs, running=False)
    running_jobs = get_running_jobs(all_jobs)
    
    completed_jobs = get_all_jobs(completed_pattern)
    remaining_jobs = {k: v for k, v in all_jobs.items()
                      if k not in completed_jobs and k not in queued_jobs}

    print('{} total jobs.'.format(len(all_jobs)))
    print('{} completed jobs.'.format(len(completed_jobs)))
    print('{} queued jobs.'.format(len(queued_jobs)))
    print('{} running jobs.'.format(len(running_jobs)))
    print('{} remaining jobs.'.format(len(remaining_jobs)))

    to_submit = min(len(remaining_jobs), max(0, n_jobs-len(queued_jobs)))
    print('Submitting {} more jobs.'.format(to_submit))

    for name, cwd in sorted(remaining_jobs.items())[:to_submit]:
        cmd = 'sbatch -p {0} -t {1} --wrap="sh {2}.sh" -J {2}'.format(queue, time, name)
        print(cwd, cmd)
        subprocess.run(cmd, cwd=cwd, shell=True)

################################################################################

template_stage1="""INPUT_FILE  {grid}_out.mae

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
"""

@main.command()
@click.argument('root')
def setup_stage1(root):
    from schrodinger.structure import StructureReader
    df = pd.read_csv('{}/structures/pdb.csv'.format(root))
    grid = [fname for fname in os.listdir('{}/docking/grids'.format(root)) if fname[0] != '.']
    assert len(grid) ==  1
    grid = grid[0]
    print(grid)
    for _, row in df.iterrows():
        print(row['ID'])
        _setup_stage1(row['ID']+'_lig', grid, root)

@main.command()
@click.argument('cwd')
def check_stage1(cwd):
    print(cwd, _check_stage1(cwd))

@main.command()
@click.argument('workdir')
def clear_stage1(workdir):
    if not _check_stage1(workdir):
        print(workdir)


def _setup_stage1(ligand, grid, root):
    from schrodinger.structure import StructureReader
    name = '{}-to-{}'.format(ligand, grid)
    cwd = '{}/docking/ifd2/{}'.format(root, name)

    cmd = ('{}/ifd -NGLIDECPU 1 -NPRIMECPU 1 {}-stage1.inp -HOST localhost'
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

    with open('{}/{}-stage1.inp'.format(cwd, name), 'w') as fp:
        fp.write(template_stage1.format(grid=grid, ligand=ligand, resid=resid))

    with open('{}/{}-stage1.sh'.format(cwd, name), 'w') as fp:
        fp.write(cmd + '\n')

def _clear_stage1(cwd, name):
    os.system('rm -rf {}/{}-stage1*'.format(cwd, name))
    os.system('rm slurm*')
    os.system('rm sh*')

def _attempted_stage1(cwd, name):
    return os.path.exists('{}/{}-stage1.log'.format(cwd, name))

def _get_stage1_workdir(cwd):
    workdir = sorted(glob(cwd + '/*-stage1_workdir/initial_docking_dir*'))
    if not workdir:
        return False
    return workdir[-1]

def _check_stage1(cwd):
    workdir = _get_stage1_workdir(cwd)
    if not workdir:
        return False

    failed = []
    for i, a in enumerate(['A','B','C','D','E','F','G']):
        pv  = wildcard('{}/*scale_lig1_G_batchglide_0000{}_pv.maegz'.format(workdir, i))
        log = wildcard('{}/*scale_lig1_{}.log'.format(workdir, a))

        if pv: continue

        if log:
            with open(log) as fp:
                txt = fp.read()

            phrases = ['** NO ACCEPTABLE LIGAND POSES WERE FOUND **',
                       'NO VALID POSES AFTER MINIMIZATION: SKIPPING.']

            if any(phrase in txt for phrase in phrases):
                continue
        failed += [a]
    return not bool(failed)

################################################################################

template_stage2="""{inputs}

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
"""

@main.command()
@click.argument('root')
def setup_stage2(root):
    for cwd in glob(root+'/*'):
        name = cwd.split('/')[-1]
        ligand = name.split('-to-')[0]

        if os.path.exists('{}/{}-stage2.inp'.format(cwd, name)):
            print('input file exists. not overwriting.')
            continue

        if not _check_stage1(cwd):
            continue
        
        stage1_dir = _get_stage1_workdir(cwd)
        assert stage1_dir
        inputs = []
        for pv in glob('{}/*-scale_lig1_*_pv-*.maegz'.format(stage1_dir)):
            inputs += ['INPUT_FILE  {}'.format(pv)]
        inputs = '\n'.join(sorted(inputs))

        with open('{}/{}-stage2.inp'.format(cwd, name), 'w') as fp:
            fp.write(template_stage2.format(inputs=inputs, ligand=ligand))

        cmd = ('{}/ifd -NGLIDECPU 1 -NPRIMECPU 1 {}-stage2.inp -HOST localhost'
               ' -SUBHOST localhost -WAIT -STRICT -RESTART'
               ).format(os.environ['SCHRODINGER'], name)

        with open('{}/{}-stage2.sh'.format(cwd, name), 'w') as fp:
            fp.write(cmd + '\n')

@main.command()
@click.argument('workdir')
def check_stage2(workdir):
    print(workdir, _check_stage2(workdir))

def _get_stage2_workdir(cwd):
    workdir = sorted(glob(cwd + '/*-stage2_workdir/glide_docking_dir*'))
    if not workdir:
        return False
    return workdir[-1]

def _check_stage2(cwd):
    workdir = _get_stage2_workdir(cwd)
    if not workdir: return False

    stage1_workdir = _get_stage1_workdir(cwd)
    assert stage1_workdir

    stage1s = [pv.split('/')[-1].split('.')[0]
               for pv in glob('{}/*pv-*.maegz'.format(stage1_workdir))]

    failed = []
    for stage1 in stage1s:
        pv  = '{}/{}-1_pv-1.maegz'.format(workdir, stage1)
        log = '{}/{}-1.log'.format(workdir, stage1)
        
        if os.path.exists(pv):
            continue

        if os.path.exists(log):
            with open(log) as fp:
                txt = fp.read()

            phrases = ['** NO ACCEPTABLE LIGAND POSES WERE FOUND **',
                       'NO VALID POSES AFTER MINIMIZATION: SKIPPING.']

            if any(phrase in txt for phrase in phrases):
                continue
        failed += [stage1]
    return not bool(failed)

def _clear_stage2(cwd):
    to_delete = ' '.join(glob('{}/*stage2*'.format(cwd)))
    print(to_delete)
    os.system('rm -rf {}'.format(to_delete))


################################################################################

template_stage3="""{inputs}
STAGE SCORING
  SCORE_NAME  r_psp_IFDScore
  TERM 1.000,r_psp_Prime_Energy,1
  TERM 9.057,r_i_glide_gscore,0
  TERM 1.428,r_i_glide_ecoul,0
  REPORT_FILE report.csv
"""

@main.command()
@click.argument('root')
def setup_stage3(root):
    for cwd in glob(root+'/*'):
        name = cwd.split('/')[-1]
        ligand = name.split('-to-')[0]

        if os.path.exists('{}/{}-stage3.inp'.format(cwd, name)):
            print('input file exists. not overwriting.')
            continue

        if not _check_stage2(cwd):
            continue
        
        stage2_dir = _get_stage2_workdir(cwd)
        assert stage2_dir
        inputs = []
        for pv in glob('{}/*-1_pv-1.maegz'.format(stage2_dir)):
            inputs += ['INPUT_FILE  {}'.format(pv)]
        inputs = '\n'.join(sorted(inputs))

        with open('{}/{}-stage3.inp'.format(cwd, name), 'w') as fp:
            fp.write(template_stage3.format(inputs=inputs))

        cmd = ('{}/ifd -NGLIDECPU 1 -NPRIMECPU 1 {}-stage3.inp -HOST localhost'
               ' -SUBHOST localhost -WAIT -STRICT -RESTART'
               ).format(os.environ['SCHRODINGER'], name)

        with open('{}/{}-stage3.sh'.format(cwd, name), 'w') as fp:
            fp.write(cmd + '\n')

################################################################################

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

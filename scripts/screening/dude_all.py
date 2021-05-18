import os
import sys
import subprocess
from glob import glob
import pandas as pd
from math import ceil

"""
for i in *; do cd $i; sbatch -p rondror --wrap="python ~/combind/scripts/screening/dude_all.py $i"; cd -; done;
for i in *; do cd $i; python ~/combind/scripts/screening/dude_all.py $i; cd -; done;

for i in *; do cd $i; sbatch -p rondror --wrap="python ~/combind/scripts/screening/dude_all.py $i score"; cd -; done;
for i in *; do cd $i; python ~/combind/scripts/screening/dude_all.py $i score; cd -; done;

for i in *; do cd $i; sbatch -p rondror --wrap="$SCHRODINGER/utilities/structcat -o all/all.maegz all/ligands/[0-9][0-9][0-9][0-9][0-9].maegz"; cd -; done;
for i in *; do cd $i; sbatch -p rondror -J sim             --wrap="python $COMBINDHOME/scripts/screening/run_dude.py similarity all.smi scores/all_rd1_shape_5"; cd -; done;
for i in *; do cd $i; sbatch -p rondror -J sha -t 12:00:00 --wrap="python $COMBINDHOME/scripts/screening/run_dude.py shape all/all.maegz scores/all_rd1_shape_5"; cd -; done;
"""

def srun(cmd):
    print(cmd)
    p = subprocess.run(cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    print(p.stdout.strip())
    print(p.stderr.strip())
    if 'QOS policy' in p.stderr:
        print('Detected sbatch submission error. Exiting.')
        exit()

def split_ligands(ligand_fname, root, max_ligands=500):
    assert os.path.exists(root)

    ligands = pd.read_csv(ligand_fname, sep=' ')
    n_files = ceil(ligands.shape[0] / max_ligands)

    for i in range(n_files):
        start = i*max_ligands
        end = (i+1)*max_ligands
        fname = '{}/{}.smi'.format(root, str(i).rjust(5, '0'))

        if not os.path.exists(fname):
            ligands[start:end].to_csv(fname, index=False, sep=' ')

def ligprep(root, protein):
    for smiles in glob(root + '/*.smi'):
        if 'dropped' in smiles: continue
        mae = smiles.replace('.smi', '.maegz')
        name = os.path.basename(smiles).split('.')[0]
        if not os.path.exists(mae):
            cmd = 'sbatch -p owners -J lp_{}_{} --wrap="combind ligprep {} {} --screen" --dependency=singleton'
            cmd = cmd.format(protein, name, smiles, root)
            print(cmd)
            subprocess.run(cmd, shell=True)

def dock(ligand_root, dock_root, protein):
    for smiles in glob(ligand_root + '/*.smi'):
        if 'dropped' in smiles: continue
        mae = smiles.replace('.smi', '.maegz')
        name = os.path.basename(smiles).split('.')[0]
        pv = '{}/{}-to-XTAL/{}-to-XTAL_pv.maegz'.format(dock_root, name, name)
        
        if os.path.exists(mae) and not os.path.exists(pv):
            cmd = 'sbatch -p owners -J dock_{}_{} --wrap="combind dock {} {} --screen" -t 12:00:00 --dependency=singleton'
            cmd = cmd.format(protein, name, dock_root, mae)
            print(cmd)
            subprocess.run(cmd, shell=True)

def featurize_single(dock_root, protein):
    for pv in glob(dock_root + '/*/*_pv.maegz'):
        ifp = pv.replace('_pv.maegz', '_ifp_rd1.csv')
        name = os.path.basename(pv).split('-to-')[0]
        if not os.path.exists(ifp):
            cmd = 'sbatch -p owners -J f_{}_{} --wrap="combind featurize . {} --max-poses 100000000" -t 12:00:00 --mem=32GB --dependency=singleton'
            cmd = cmd.format(protein, name, pv)
            print(cmd)
            subprocess.run(cmd, shell=True)

def featurize(protein, scores):
    cwd = os.getcwd()
    os.chdir(scores)
    for pv in glob('../../../all/docking/*/*_pv.maegz'):
        name = os.path.basename(pv).split('-to-')[0]
        ifp = pv.replace('_pv.maegz', '_ifp_rd1.csv')
        shape = 'screen/shape/shape-{}-to-XTAL_pv-and-binder_pv.npy'.format(name)
        ifp_pair = ['screen/ifp-pair/{}-{}-to-XTAL_ifp_rd1-and-binder_ifp_rd1.npy'.format(feature, name)
                    for feature in ['hbond', 'saltbridge', 'contact']]
        done = all(os.path.exists(x) for x in ifp_pair)
        done = done and os.path.exists(shape)
        if os.path.exists(ifp) and not done:
            cmd = ('combind featurize screen '
                   '{} bpp/binder_pv.maegz --no-mcss --screen')
            cmd = cmd.format(pv)
            cmd = 'sbatch -p owners -t 12:00:00 -J f2_{}_{} --dependency=singleton --wrap="{}"'.format(protein, name, cmd)
            srun(cmd)
    os.chdir(cwd)

def verify_featurize(protein, scores):
    cwd = os.getcwd()
    os.chdir(scores)
    for pv in glob('../../../all/docking/*/*_pv.maegz'):
        name = os.path.basename(pv).split('-to-')[0]
        anno_pv = 'screen/anno/{}_pv.maegz'.format(os.path.basename(pv).split('_pv')[0])
        if not os.path.exists(anno_pv):
            shape = 'screen/shape/shape-{}-to-XTAL_pv-and-binder_pv.npy'.format(name)
            ifp = pv.replace('_pv.maegz', '_ifp_rd1.csv')
            cmd = ('combind featurize screen --verify --delete '
                   '{} bpp/binder_pv.maegz --no-mcss --screen')
            cmd = cmd.format(pv)
            srun(cmd)
    os.chdir(cwd)

def screen(protein, features, scores, stats='/oak/stanford/groups/rondror/users/jpaggi/dude_stats'):
    cwd = os.getcwd()
    os.chdir(scores)
    for name in glob('screen/shape/shape*.npy'):
        name = os.path.basename(name).split('_pv-and-')[0].replace('shape-', '')
        
        if not os.path.exists('screen/scores'):
            os.mkdir('screen/scores')
        if not os.path.exists('screen/anno'):
            os.mkdir('screen/anno')
        if not os.path.exists('screen/combind'):
            os.mkdir('screen/combind')
        if not os.path.exists('screen/glide'):
            os.mkdir('screen/glide')

        score_fname = 'screen/scores/{}.npy'.format(name)
        orig_pv = '../../../all/docking/{}/{}_pv.maegz'.format(name, name)
        anno_pv = 'screen/anno/{}_pv.maegz'.format(name)
        combind_top = 'screen/combind/{}_pv.maegz'.format(name)
        glide_top = 'screen/glide/{}_pv.maegz'.format(name)

        if not os.path.exists(score_fname):
            cmd = ("combind screen {} "
                   "../../../all/docking/{}/{}_gscore.npy "
                   "--ifp-fname screen/ifp-pair/{}-{}_ifp_rd1-and-binder_ifp_rd1.npy "
                   "--shape-fname screen/shape/shape-{}_pv-and-binder_pv.npy "
                   "--stats-root {}/{} "
                   "--features {}")
            cmd = cmd.format(score_fname, name, name, '{}', name, name, stats, protein, features)
            print(cmd)
            subprocess.run(cmd, shell=True)

        if os.path.exists(score_fname) and not os.path.exists(anno_pv):
            cmd = 'combind apply-scores {} {} {}'
            cmd = cmd.format(orig_pv, score_fname, anno_pv)
            print(cmd)
            subprocess.run(cmd, shell=True)

        if os.path.exists(anno_pv) and not os.path.exists(combind_top):
            cmd = '$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score -o {} {}'
            cmd = cmd.format(combind_top, anno_pv)
            print(cmd)
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')

            if 'Incomplete CT block' in proc.stderr:
                print('Removing {}'.format(anno_pv))
                os.remove(anno_pv)

        if os.path.exists(anno_pv) and not os.path.exists(glide_top):
            cmd = '$SCHRODINGER/utilities/glide_sort -best_by_title -o {} {}'
            cmd = cmd.format(glide_top, anno_pv)
            print(cmd)
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')

            assert 'Incomplete CT block' not in proc.stderr

    os.chdir(cwd)

def merge(scores):
    cwd = os.getcwd()
    os.chdir(scores)

    names = glob('../../../all/ligands/[0-9][0-9][0-9][0-9][0-9].smi')
    names = [os.path.basename(name).split('.')[0] for name in names]

    combind_pv = 'screen/combind_pv.maegz'
    combind_csv = 'combind.csv'
    glide_pv = 'screen/glide_pv.maegz'
    glide_csv = 'glide.csv'

    if not os.path.exists(combind_pv):
        fnames = ['screen/combind/{}-to-XTAL_pv.maegz'.format(name) for name in names]
        ready = all([os.path.exists(fname) for fname in fnames])
        if ready:
            cmd = '$SCHRODINGER/utilities/glide_sort -best_by_title -use_prop_d r_i_combind_score -o {} {}'
            cmd = cmd.format(combind_pv, ' '.join(fnames))
            print(cmd)
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True, encoding='utf-8')

            for line in proc.stderr.split('\n'):
                if 'Incomplete CT block' in line:
                    file = line.split("'")[1]
                    print('Removing {}'.format(file))
                    os.remove(file)

    if not os.path.exists(glide_pv):
        fnames = ['screen/glide/{}-to-XTAL_pv.maegz'.format(name) for name in names]
        ready = all([os.path.exists(fname) for fname in fnames])
        if ready:
            cmd = '$SCHRODINGER/utilities/glide_sort -best_by_title -o {} {}'
            cmd = cmd.format(glide_pv, ' '.join(fnames))
            print(cmd)
            subprocess.run(cmd, shell=True)

    if os.path.exists(combind_pv) and not os.path.exists(combind_csv):
        cmd = 'combind scores-to-csv {} {}'
        cmd = cmd.format(combind_pv, combind_csv)
        print(cmd)
        subprocess.run(cmd, shell=True)
    
    if os.path.exists(glide_pv) and not os.path.exists(glide_csv):
        cmd = 'combind scores-to-csv {} {}'
        cmd = cmd.format(glide_pv, glide_csv)
        print(cmd)
        subprocess.run(cmd, shell=True)

    os.chdir(cwd)

if not os.path.exists('all'):
    os.mkdir('all')
if not os.path.exists('all/ligands'):
    os.mkdir('all/ligands')
if not os.path.exists('all/docking'):
    os.mkdir('all/docking')

if len(sys.argv) == 2:
    # split_ligands('all.smi', 'all/ligands')
    # ligprep('all/ligands', sys.argv[1])
    # dock('all/ligands', 'all/docking', sys.argv[1])
    #featurize_single('all/docking', sys.argv[1])
    featurize(sys.argv[1], 'scores/all_rd1_shape_0/0')
    for n in [1, 3, 5, 10]:
        for i in range(5):
            scores = 'scores/all_rd1_shape_{}/{}'.format(n, i)
            featurize(sys.argv[1], scores)

elif len(sys.argv) == 3 and sys.argv[2] == 'verify':
    verify_featurize(sys.argv[1], 'scores/all_rd1_shape_0/0')
    for n in [1, 3, 5, 10]:
        for i in range(5):
            scores = 'scores/all_rd1_shape_{}/{}'.format(n, i)
            verify_featurize(sys.argv[1], scores)

elif len(sys.argv) == 3 and sys.argv[2] == 'score':
    scores = 'scores/all_rd1_shape_0/0'
    screen(sys.argv[1], 'shape,contact,hbond,saltbridge', scores)
    merge(scores)
    for n in [1, 3, 5, 10]:
        for i in range(5):
            scores = 'scores/all_rd1_shape_{}/{}'.format(n, i)
            screen(sys.argv[1], 'shape,contact,hbond,saltbridge', scores)
            merge(scores)

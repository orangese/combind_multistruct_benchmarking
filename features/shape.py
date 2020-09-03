import subprocess
import os
import tempfile
import numpy as np
from schrodinger.structure import StructureReader, StructureWriter

CMD = '$SCHRODINGER/shape_screen -shape {poses1} -screen {poses2} -inplace {typing} {norm} -distinct -NOJOBID'

def write_and_count(pv_in, pv_out, max_poses):
    n = 0
    with StructureReader(pv_in) as sts, StructureWriter(pv_out) as writer:
        next(sts)
        for i, st in enumerate(sts):
            writer.append(st)
            n += 1
            if n == max_poses:
                break
    return n

def shape(pv1, pv2, version='pharm_max', max_poses=float('inf')):
    typing, norm = version.split('_')
    
    if typing == 'pharm':
        typing = '-pharm'
    elif typing == 'mmod':
        typing = '-atomtypes mmod'
    elif typing == 'element':
        typing = '-atomtypes element'
    elif typing == 'qsar':
        typing = '-atomtypes qsar'
    else:
        assert False, 'Typing {} not supported.'.format(typing)

    if norm == 'max':
        norm = '-norm 1'
    elif norm == 'min':
        norm = '-norm 2'
    else:
        assert False, 'Norm {} not supported.'.format(norm)

    with tempfile.TemporaryDirectory() as wd:
        poses1 = wd+'/poses1.maegz'
        poses2 = wd+'/poses2.maegz'
        output = wd+'/poses1_align.maegz'

        n_poses1 = write_and_count(pv1, poses1, max_poses)
        n_poses2 = write_and_count(pv2, poses2, max_poses)

        cmd = CMD.format(poses1=os.path.basename(poses1),
                         poses2=os.path.basename(poses2),
                         typing=typing, norm=norm)
        subprocess.run(cmd, shell=True, cwd=wd)

        sims = np.zeros((n_poses1, n_poses2))
        with StructureReader(output) as sts:
            for k, st in enumerate(sts):
                i = k % n_poses1
                j = int(k / n_poses1)
                sims[i, j] = st.property['r_phase_Shape_Sim']
            assert k == n_poses1*n_poses2-1
    return sims

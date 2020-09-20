import subprocess
import os
import tempfile
import numpy as np
from schrodinger.structure import StructureReader, StructureWriter

CMD = '$SCHRODINGER/shape_screen -shape {poses1} -screen {poses2} -inplace {typing} {norm} -distinct -NOJOBID'


def key(st):
    if 'i_i_glide_posenum' not in st.property:
        return st.title
    return (st.property['i_m_source_file_index'],
            st.property['i_i_glide_posenum'],
            st.property['s_lp_Variant'])


def write_and_count(pv_in, pv_out, max_poses):
    titles = []
    with StructureReader(pv_in) as sts, StructureWriter(pv_out) as writer:
        next(sts)
        for i, st in enumerate(sts):
            writer.append(st)
            titles += [key(st)]
            if len(titles) == max_poses:
                print('max reached')
                break
    assert len(set(titles)) == len(titles)
    return titles

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
        log    = wd+'/poses1_shape.log'

        ligands1 = write_and_count(pv1, poses1, max_poses)
        ligands2 = write_and_count(pv2, poses2, max_poses)


        cmd = CMD.format(poses1=os.path.basename(poses1),
                         poses2=os.path.basename(poses2),
                         typing=typing, norm=norm)
        subprocess.run(cmd, shell=True, cwd=wd)

        if not os.path.exists(output):
            with open(log) as fp:
                txt = fp.read()
            assert 'Reference shape must contain at least 3 spheres' in txt, txt
            return 0.5*np.ones((len(ligands1), len(ligands2)))

        sims = np.zeros((len(ligands1), len(ligands2)))
        with StructureReader(output) as sts:
            for k, st in enumerate(sts):
                i = k % len(ligands1)
                j = ligands2.index(key(st))
                sims[i, j] = st.property['r_phase_Shape_Sim']
    return sims

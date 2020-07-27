import subprocess
import os
import tempfile
import numpy as np
from schrodinger.structure import StructureReader, StructureWriter

CMD = '$SCHRODINGER/shape_screen -shape {pose1} -screen {poses2} -inplace {typing} {norm} -distinct -NOJOBID'

def shape(pv1, pv2, version='pharm_max', max_poses=100):
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
        pose1  = wd+'/{}.maegz'
        poses2 = wd+'/poses2.maegz'
        output = wd+'/{}_align.maegz'

        # write first pv pose by pose to files
        n_poses1 = 0
        with StructureReader(pv1) as sts:
            next(sts)
            for i, st in enumerate(sts):
                st.write(pose1.format(i))
                n_poses1 += 1
                if n_poses1 == max_poses:
                    break

        # write out pose only version of pv2
        n_poses2 = 0
        with StructureReader(pv2) as sts, StructureWriter(poses2) as writer:
            next(sts)
            for st in sts:
                writer.append(st)
                n_poses2 += 1
                if n_poses2 == max_poses:
                    break

        # execute shape_sim
        for i in range(n_poses1):
            cmd = CMD.format(pose1=os.path.basename(pose1.format(i)),
                             poses2=os.path.basename(poses2),
                             typing=typing, norm=norm)
            subprocess.run(cmd, shell=True, cwd=wd)

        # read results
        sims = np.zeros((n_poses1, n_poses2))
        for i in range(n_poses1):
            with StructureReader(output.format(i)) as sts:
                for j, st in enumerate(sts):
                    sims[i, j] = st.property['r_phase_Shape_Sim']
    return sims

if __name__ == '__main__':
    import click
    @click.command()
    @click.argument('pv1')
    @click.argument('pv2')
    @click.argument('output')
    @click.option('--version', default='pharm_max')
    @click.option('--max-poses', default=100)
    def main(pv1, pv2, output, version, max_poses):
        assert 'npy' in output, 'Output is in npy format.'
        sims = shape(pv1, pv2, version=version, max_poses=max_poses)
        np.save(output, sims)
    main()

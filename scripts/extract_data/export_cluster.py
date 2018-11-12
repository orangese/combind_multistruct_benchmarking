import os
import sys
from schrodinger.structure import StructureReader, StructureWriter

#print sys.argv

glide_dir = sys.argv[2] # '/scratch/PI/rondror/CREB_inhibitors/glide_output3'
output_dir = sys.argv[3] # '/scratch/PI/rondror/CREB_inhibitors/analysis/output_clusters'

def export(poses_to_export, cluster_title):
    os.system('mkdir -p {}'.format(output_dir))
    if os.path.exists('{}/{}.mae'.format(output_dir, cluster_title)):
        os.system('rm {}/{}.mae'.format(output_dir, cluster_title))
    cluster_out = StructureWriter('{}/{}.mae'.format(output_dir, cluster_title))

    for j, (docked_pair, pnum) in enumerate(poses_to_export):
        #title = '{}-{}'.format(pnum, docked_pair.split('_')[0])
        print docked_pair, pnum 
        if pnum == 'T': 
        #path = '{}/{}/score_pv.maegz'.format(glide_dir, docked_pair)
            path = glide_dir.split('/')
            path[-1] = 'ligands'
            path[-2] = 'structures'
            path.append('{}_lig.mae'.format(docked_pair))
            path = '/'.join(path)
        elif pnum == 'R': path = '{}/{}/refine_pv.maegz'.format(glide_dir, docked_pair)
        else: path = '{}/{}/{}_pv.maegz'.format(glide_dir, docked_pair, docked_pair)

        if pnum in ['T','R']: st_n = 0
        else: st_n = pnum
        
        for i, st in enumerate(StructureReader(path)):
            if j == 0 and i == 0: pass
                #st._setTitle('prot-{}'.format(docked_pair.split('-to-')[1]))
                #cluster_out.append(st)
            if i - 1 == st_n:
                st._setTitle('{}-{}'.format(pnum, docked_pair.split('_')[0]))
                cluster_out.append(st)
                break

        #if pnum == 'T':
        #    for i, st in enumerate(StructureReader('{}/{}/score_pv.maegz'.format(glide_dir, docked_pair)):
        #        if i == 1:
        #            pose = st
        #elif pnum == 'R':
        #    for i, st in enumerate(StructureReader('{}/{}/refine_pv.maegz'.format(glide_dir, docked_pair)):
        #        if i == 1:
        #            pose = st
        #else:
        #    for i, st in enumerate(StructureReader('{}/{}/{}_pv.maegz'.format(glide_dir, docked_pair, docked_pair))):
        #        if i - 1 == pnum:
        #            st._setTitle('{}-{}'.format(pnum, docked_pair))
        #            cluster_out.append(st)
        #            break
        #pose._setTitle(title)
        #cluster_out.append(pose)
    cluster_out.close()

cluster_name = sys.argv[1]

cluster = []
for i, item in enumerate(sys.argv[4:]):
    if i % 2 == 0:
        pose = sys.argv[4:][i+1]
        if pose in ['T','R']: cluster.append((item, pose))
        else: cluster.append((item, int(pose)))
export(cluster, cluster_name)


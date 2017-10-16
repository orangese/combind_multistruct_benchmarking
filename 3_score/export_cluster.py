import os
import sys
from schrodinger.structure import StructureReader, StructureWriter

glide_dir = sys.argv[2] # '/scratch/PI/rondror/CREB_inhibitors/glide_output3'
output_dir = sys.argv[3] # '/scratch/PI/rondror/CREB_inhibitors/analysis/output_clusters'

def export(name_to_pnum, cluster_title):
    os.system('mkdir -p {}'.format(output_dir))
    if os.path.exists('{}/{}.mae'.format(output_dir, cluster_title)):
        os.system('rm {}/{}.mae'.format(output_dir, cluster_title))
    cluster_out = StructureWriter('{}/{}.mae'.format(output_dir, cluster_title))
    for docked_pair, pnum in name_to_pnum.items():
        for i, st in enumerate(StructureReader('{}/{}/{}_pv.maegz'.format(glide_dir, docked_pair, docked_pair))):
            if i - 1 == pnum:
                st._setTitle('{}-{}'.format(pnum, docked_pair))
                cluster_out.append(st)
                break
    cluster_out.close()

cluster_name = sys.argv[1]

cluster = {}
for i, item in enumerate(sys.argv[4:]):
    if i % 2 == 0:
        cluster[item] = int(sys.argv[4:][i+1])

export(cluster, cluster_name)


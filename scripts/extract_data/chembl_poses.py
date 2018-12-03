import os
import sys
sys.path.append('../../dock')
sys.path.append('../../2_ifp')
sys.path.append('../../score')

from containers import Dataset
from shared_paths import shared_paths
from schrodinger.structure import StructureReader, StructureWriter

def read_score_file(fname):
    """
    Returns a dictionary of {lig_name: pose, ...} stored
    in FNAME.
    """
    pose_cluster = {}
    with open(fname) as fp:
        for line in fp:
            tok = line.strip().split(',')
            if tok[0] != 'max_score':
                pose_cluster[tok[0]] = int(tok[1])
    return pose_cluster

scores_version, protein, ligand = sys.argv[1:]
out = '/scratch/PI/rondror/combind/bpp_outputs/chembl/{}_{}_{}'.format(scores_version.replace('/', '-'), protein, ligand)

data = Dataset(shared_paths, [protein])
struct = data.proteins[protein].lm.st
docking =  data.proteins[protein].docking[struct]

scores_path = "{}/{}/scores/{}/{}-to-{}.sc".format(shared_paths['data'], protein, scores_version, ligand, struct)
pose_cluster = read_score_file(scores_path)


data.load({protein: list(set([lig for lig in pose_cluster]))}, load_mcss = False, load_fp = False)

pv_template = "{0:}/{1:}/docking/{2:}/{3:}-to-{4:}/{3:}-to-{4:}_pv.maegz".format(shared_paths['data'],
                                                                                 protein,
                                                                                 shared_paths['docking'],
                                                                                 '{0:}', struct)
            
with StructureWriter(out+'_glide_pv.mae') as glide, \
     StructureWriter(out+'_combind_pv.mae') as combind:
    for i, (ligand, pose) in enumerate(pose_cluster.items()):
        poses = docking.ligands[ligand].poses
        pv = list(StructureReader(pv_template.format(ligand)))
        grid = pv[0]
        if not i:
            glide.append(grid)
            combind.append(grid)
        glide.append(pv[1])
        combind.append(pv[pose+1])

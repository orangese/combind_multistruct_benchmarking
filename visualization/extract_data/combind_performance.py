"""
Ideal behaviour would be to
- Tabulate performance for all datasets
- Create poseviewer files for top glide + combind poses
- Get fingerprints for all top poses?
- Get energy terms?
"""
import os
import sys
from containers import Dataset
sys.path.append('../1_dock')
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

if __name__ == '__main__':
    scores_version = sys.argv[1]
    base = '/scratch/PI/rondror/combind/bpp_outputs/'+scores_version
    with open('{}/rmsds.csv'.format(base), 'w') as fp:
        for protein in os.listdir(shared_paths['data']):
            print(protein)
            if protein[0] == '.': continue
            scores_path = "{}/{}/scores/{}/".format(shared_paths['data'], protein, scores_version)
            if not os.path.exists(scores_path): continue
            pose_clusters = [read_score_file(scores_path+fname)
                             for fname in os.listdir(scores_path) if fname[-3:] == '.sc']
            data = Dataset(shared_paths, [protein])
            data.load({protein: list(set([lig for pose_cluster in pose_clusters for lig in pose_cluster]))}, load_mcss = False, load_fp = False)
            
            # print RMSDs
            struct = data.proteins[protein].lm.st
            docking =  data.proteins[protein].docking[struct]
            
            for pose_cluster in pose_clusters:
                for lig, pose in pose_cluster.items():
                    if 'CHEMBL' in lig: continue
                    poses = docking.ligands[lig].poses
                    best_rmsd = min(map(lambda x: x.rmsd, poses[:100]))
                    combind_rmsd = poses[pose].rmsd
                    glide_rmsd = poses[0].rmsd
                    fp.write(','.join(map(str, [protein, lig, best_rmsd, glide_rmsd, combind_rmsd]))+'\n')

            # Write poseviewers of top poses
            pv_template = "{0:}/{1:}/docking/{2:}/{3:}-to-{4:}/{3:}-to-{4:}_pv.maegz".format(shared_paths['data'],
                                                                                             protein,
                                                                                             shared_paths['docking'],
                                                                                             '{0:}', struct)
            with StructureWriter('{}/glide_{}_pv.mae'.format(base, protein)) as glide, \
                 StructureWriter('{}/combind_{}_pv.mae'.format(base, protein)) as combind, \
                 StructureWriter('{}/best_{}_pv.mae'.format(base, protein)) as best, \
                 StructureWriter('{}/crystal_{}_pv.mae'.format(base, protein)) as crystal:
        
                for i, pose_cluster in enumerate(pose_clusters):
                    for ligand, pose in pose_cluster.items():
                        if 'CHEMBL' in ligand: continue
                        poses = docking.ligands[ligand].poses
                        best_pose, best_rmsd = None, float('inf')
                        for p in poses:
                            if p.rmsd < best_rmsd:
                                best_pose, best_rmsd = p.rank, p.rmsd
                        pv = list(StructureReader(pv_template.format(ligand)))
                        grid = pv[0]
                        if not i:
                            glide.append(grid)
                            combind.append(grid)
                            best.append(grid)
                        glide.append(pv[1])
                        combind.append(pv[pose+1])
                        best.append(pv[best_pose+1])
                        if ligand == '1E1V_lig': continue
                        crystal.extend(StructureReader("{}/{}/structures/ligands/{}.mae".format(shared_paths['data'], protein, ligand)))

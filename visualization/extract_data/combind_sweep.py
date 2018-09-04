import os
import sys
sys.path.append('../../1_dock')
sys.path.append('../../2_ifp')
sys.path.append('../../3_analyze')

from containers import Dataset
from shared_paths import shared_paths

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
            data = Dataset(shared_paths, [protein])
            for settings in os.listdir(scores_path):
                pose_clusters = [read_score_file(scores_path+settings+'/'+fname)
                                 for fname in os.listdir(scores_path+settings) if fname[-3:] == '.sc']
                

            
                # print RMSDs
                struct = data.proteins[protein].lm.st
                docking =  data.proteins[protein].docking[struct]
            
                for pose_cluster in pose_clusters:
                    for lig, pose in pose_cluster.items():
                        if 'CHEMBL' in lig: continue
                        if lig  not in docking.ligands:
                            data.load({protein: [lig]}, load_mcss = False, load_fp = False)
                        poses = docking.ligands[lig].poses
                        best_rmsd = min(map(lambda x: x.rmsd, poses[:100]))
                        combind_rmsd = poses[pose].rmsd
                        glide_rmsd = poses[0].rmsd
                        fp.write(','.join(map(str, [settings, protein, lig, best_rmsd, glide_rmsd, combind_rmsd]))+'\n')

import os
import sys

from containers import Protein
from statistics import statistics
from prob_opt import PredictStructs

class ScoreContainer:
    """
    Manages execution of ComBind for a given protein and settings.
    """
    def __init__(self, root, prot, struct):
        self.prot = prot
        self.struct = struct
        self.root = root
      
        self.settings = self.read_settings()
        self.sp = self.settings['shared_paths']
 
        self.stats = self.init_stats()
        self.predict_data = Protein(prot)
        self.ps = PredictStructs(self.predict_data, self.stats, 
                                 self.settings['k_list'], self.settings['num_poses'],
                                 self.settings['t'])

    def read_settings(self):
        tr = {}
        with open('{}/settings.py'.format(self.root)) as f:
            for line in f:
                var,val = line.split('=')
                tr[var] = eval(val)
        return tr

    def init_stats(self):
        data = {}
        for protein in self.settings['stats_prots']:
            lm = Protein(protein).lm
            ligands = lm.docked(lm.pdb)[:self.settings['num_stats_ligs']+1]
            self_docked = lm.st+'_lig'
            if self_docked in ligands:
                ligands.remove(self_docked)
            else:
                ligands.pop(-1)
            data[protein] = ligands
        return statistics(data, self.settings['k_list'])

    def compute_results_chembl(self, query):
        assert self.settings['chembl']
        chembl_ligs = self.predict_data.lm.get_helpers(query, self.settings['chembl_file'],
                                          num=self.settings['num_pred_chembl'],
                                          struct=self.struct)
        return self.compute_results([query]+chembl_ligs)

    def compute_results(self, queries):
        self.predict_data.load_docking(queries, self.struct)
        best_cluster, all_scores, all_rmsds = self.ps.max_posterior(queries, restart=15,
                                                                    sampling=3)
        return best_cluster

    def write_results(self, cluster, fname):
        """
        Part 1:
        ligand id, combind rank, combind rmsd, glide rank, glide rmsd, best rank, best rmsd
        Part 2:
        combind score, glide score, best score (if all pdb, else 0)
        """
        with open(fpath, 'w') as f:
            f.write('lig,combind_rank,combind_rmsd,glide_rank,glide_rmsd,best_rank,best_rmsd\n')
            best_cluster = {}
            for lig, combind_pose in sorted(cluster.items()):
                poses = self.predict_data.proteins[self.prot].docking[self.struct].ligands[lig].poses
                best_rmsd = float('inf')
                for i, pose in enumerate(poses[:self.settings['num_poses']]):
                    if pose.rmsd is not None and pose.rmsd < best_rmsd:
                        best_cluster[lig] = i
                        best_rmsd = pose.rmsd
                f.write(','.join(map(str, [lig,
                                           combind_pose, poses[combind_pose].rmsd,
                                           0, poses[0].rmsd,
                                           best_cluster[lig] if lig in best_cluster else None, best_rmsd]))+'\n')
            f.write('combind={},glide={},best={}\n'.format(
                                        sc.ps.log_posterior(cluster),
                                        sc.ps.log_posterior({k:0 for k in cluster.keys()}),
                                        sc.ps.log_posterior(best_cluster) if len(best_cluster) == len(cluster) else 0))

if __name__ == '__main__':
    struct, protein = sys.argv[1:3]
    queries = sys.argv[3:]

    sc = ScoreContainer(os.getcwd(), protein, struct)

    # Write out stats
    for dist, interactions in sc.stats.items():
        for interaction, de in interactions.items():
            with open('{}/{}_{}.txt'.format(sc.root, dist, interaction), 'w') as fp:
                fp.write(str(de)+'\n')

    # Compute results
    if sc.settings['chembl']:
        for query in queries:
            combind_cluster = sc.compute_results_chembl(query)
            fpath = '{}/{}-to-{}.sc'.format(sc.root, query, struct)
            sc.write_results(combind_cluster, fpath)
    else:
        combind_cluster = sc.compute_results(queries)
        fpath = '{}/pdb.sc'.format(sc.root)
        sc.write_results(combind_cluster, fpath)

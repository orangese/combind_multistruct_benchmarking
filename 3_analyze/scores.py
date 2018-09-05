import os
import sys

from containers import Dataset, LigandManager
from statistics import Statistics, readf
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
        self.predict_data = Dataset(self.sp, [prot], {prot:struct})
        self.ps = PredictStructs(self.predict_data.proteins[prot], self.stats, 
            self.settings['k_list'], self.settings['num_poses'], self.settings['t'], self.settings['score_mode'])

    def read_settings(self):
        tr = {}
        with open('{}/settings.py'.format(self.root)) as f:
            for line in f:
                var,val = line.split('=')
                tr[var] = eval(val)
        return tr

    def init_stats(self):
        if os.path.exists('{}/stats_{}.txt'.format(self.root, self.settings['k_list'][0])):
            stats = Statistics({},{},self.settings['k_list'])
            for k in self.settings['k_list']:
                tmp = readf('{}/stats_{}.txt'.format(self.root, k), k, stats.ind)
                for i in stats.ind:
                    stats.dist[i][k] = tmp[i]
        else:
            stats_ligs = {}
            stats_st = {}
            for p in self.settings['stats_prots']:
                lm = LigandManager(self.sp,p)
                stats_ligs[p] = lm.docked(lm.pdb)[:self.settings['num_stats_ligs']]
                stats_st[p] = lm.st
            stats = Statistics(stats_ligs, stats_st, self.settings['k_list'])
            stats.read(self.sp['data'], self.sp['stats'])
        return stats

    def compute_results_chembl(self, query):
        assert self.settings['chembl']
        prot = self.predict_data.proteins[self.prot]
        chembl_ligs = prot.lm.get_similar(query, self.settings['chembl_file'], 
                                          num=self.settings['num_pred_chembl'],
                                          struct=self.struct)
        return self.compute_results([query]+chembl_ligs)

    def compute_results(self, queries):
        self.predict_data.load({self.prot: queries},{self.prot: [self.struct]})
        best_cluster, all_scores, all_rmsds = self.ps.max_posterior(queries, restart=15, sampling=3)
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
                for pose in poses[:self.settings['num_poses']]:
                    if pose.rmsd is not None and pose.rmsd < best_rmsd:
                        best_cluster[lig] = pose.rank
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
    print('exit')

    # Write out stats
    for k in sc.stats.k_list:
        out_f = '{}/stats_{}.txt'.format(sc.root, k)
        sc.stats.write(out_f, k)

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

import os
import sys

from containers import Protein
from statistics import statistics
from prob_opt import PredictStructs
from density_estimate import DensityEstimate

class ScoreContainer:
    """
    Manages execution of ComBind for a given protein and settings.
    """
    def __init__(self, root, stats_root, prot, struct):
        self.prot = prot
        self.struct = struct
        self.root = root
      
        self.settings = self.read_settings()
        self.sp = self.settings['shared_paths']

        self.read_stats(stats_root)
 
        self.predict_data = Protein(prot)
        self.ps = PredictStructs(self.predict_data, self.stats,
                                 self.settings['k_list'], self.settings['num_poses'],
                                 self.settings['alpha'])

    def read_settings(self):
        tr = {}
        with open('{}/settings.py'.format(self.root)) as f:
            for line in f:
                var,val = line.split('=')
                tr[var] = eval(val)
        return tr

    def read_stats(self, stats_root):
        self.stats = {'native':{}, 'reference':{}}
        for dist in self.stats:
            for interaction in self.settings['k_list']:
                fname = '{}/{}_{}.txt'.format(stats_root, dist, interaction)
                assert os.path.exists(fname), fname
                self.stats[dist][interaction] = DensityEstimate.read(fname)

    def compute_results_chembl(self, query):
        assert self.settings['chembl']
        chembl_ligs = self.predict_data.lm.get_helpers(query, self.settings['chembl_file'],
                                          num=self.settings['num_pred_chembl'],
                                          struct=self.struct)
        return self.compute_results([query]+chembl_ligs)

    def compute_results(self, queries):
        self.predict_data.load_docking(queries, load_fp = True,
                                       load_mcss = 'mcss' in self.settings['k_list'],
                                       st = self.struct)

        if 'use_crystal_pose' in self.settings and self.settings['use_crystal_pose']:
            crystal_lig = '{}_crystal_lig'.format(self.predict_data.lm.st)
            if crystal_lig not in queries: queries += [crystal_lig]
            self.predict_data.load_docking([crystal_lig], load_crystal = True,
                                           load_fp = True,
                                           load_mcss = 'mcss' in self.settings['k_list'],
                                           st = self.struct)

        best_cluster = self.ps.max_posterior(queries, restart=50, sampling=3)
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
                poses = self.predict_data.docking[self.struct].ligands[lig].poses
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
    stats_root, struct, protein = sys.argv[1:4]
    queries = sys.argv[4:]
    sc = ScoreContainer(os.getcwd(), stats_root, protein, struct)

    if sc.settings['chembl']:
        for query in queries:
            combind_cluster = sc.compute_results_chembl(query)
            fpath = '{}/{}-to-{}.sc'.format(sc.root, query, struct)
            sc.write_results(combind_cluster, fpath)
    else:
        # This needs to be before compute results as therein queries is mutated
        fpath = '{}/{}.sc'.format(sc.root, 'pdb' if len(queries) > 1 else queries[0])
        combind_cluster = sc.compute_results(queries)
        sc.write_results(combind_cluster, fpath)

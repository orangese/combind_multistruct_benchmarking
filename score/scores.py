import os
import sys

from containers import Protein
from score.statistics import statistics
from score.prob_opt import PredictStructs
from score.density_estimate import DensityEstimate

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
        self.ps = PredictStructs({}, self.predict_data.lm.mcss, self.stats,
                                 self.settings['k_list'], self.settings['num_poses'],
                                 self.settings['alpha'])

    def read_settings(self):
        tr = {}
        with open('{}/settings.py'.format(self.root)) as f:
            for line in f:
                var, val = line.split('=')
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

        # Set ligands and optimize!
        self.ps.ligands = {lig: self.predict_data.docking[self.struct].ligands[lig]
                           for lig in queries}
        best_cluster = self.ps.max_posterior()
        return best_cluster

    def write_results(self, cluster, fname):
        """
        Part 1:
        ligand id, combind rank, combind rmsd, glide rank, glide rmsd, best rank, best rmsd
        Part 2:
        combind score, glide score, best score (if all pdb, else 0)
        """
        with open(fname, 'w') as f:
            f.write('lig,combind_rank,combind_rmsd,glide_rank,glide_rmsd,best_rank,best_rmsd\n')
            glide_cluster, best_cluster = {}, {}
            for lig, combind_pose in sorted(cluster.items()):
                poses = self.predict_data.docking[self.struct].ligands[lig].poses
                best_rmsd = float('inf')
                best_emodel = float('inf')
                for i, pose in enumerate(poses[:self.settings['num_poses']]):
                    if pose.rmsd is not None and pose.rmsd < best_rmsd:
                        best_cluster[lig] = i
                        best_rmsd = pose.rmsd
                    if pose.emodel < best_emodel:
                        glide_cluster[lig] = i
                        best_emodel = pose.emodel

                glide_pose = 0 #glide_cluster[lig]
                best_pose = best_cluster[lig] if lig in best_cluster else None
                f.write(','.join(map(str, [lig,
                                           combind_pose, poses[combind_pose].rmsd,
                                           glide_pose, poses[glide_pose].rmsd,
                                           best_pose, best_rmsd
                                           ]))+'\n')
            f.write('combind={},glide={},best={}\n'.format(
                    self.ps.log_posterior(cluster),
                    self.ps.log_posterior(glide_cluster),
                    self.ps.log_posterior(best_cluster) if len(best_cluster) == len(cluster) else 0))

    def read_results(self, fname):
        cluster = {}
        with open('{}/{}'.format(self.root, fname)) as fp:
            fp.readline()
            for line in fp:
                if line[:7] == 'combind': continue

                lig, combind_pose, *rest = line.split(',')
                cluster[lig] = int(combind_pose)

        self.predict_data.load_docking(list(cluster.keys()), load_fp = True,
                                       load_mcss = 'mcss' in self.settings['k_list'],
                                       st = self.struct)
        self.ps.ligands = {lig: self.predict_data.docking[self.struct].ligands[lig]
                           for lig in list(cluster.keys())}
        return cluster


def main(args):
    stats_root, struct, protein = args[1:4]
    queries = args[4:]
    sc = ScoreContainer(os.getcwd(), stats_root, protein, struct)

    if sc.settings['chembl']:
        for query in queries:
            combind_cluster = sc.compute_results_chembl(query)
            fname = '{}/{}-to-{}.sc'.format(sc.root, query, struct)
            sc.write_results(combind_cluster, fname)
    else:
        # This needs to be before compute results as therein queries is mutated
        fname = '{}/{}.sc'.format(sc.root, 'pdb' if len(queries) > 1 else queries[0])
        combind_cluster = sc.compute_results(queries)
        sc.write_results(combind_cluster, fname)

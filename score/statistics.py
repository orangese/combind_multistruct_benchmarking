import numpy as np
from containers import Protein
from score.density_estimate import DensityEstimate
from score.pairs import LigPair
from glob import glob
import os

class Statistics:
    def __init__(self, proteins, interactions, settings, paths, path=None):
        self.proteins = proteins
        self.interactions = interactions
        self.settings = settings
        self.distributions = ['native', 'reference']
        self.paths = paths
        self.path = path

        if self.path:
            os.makedirs(os.path.dirname(path), exist_ok=True)

        self.stats = self._load()

    @classmethod
    def read_merged(cls, path):
        interactions = set()
        stats = {}
        for fname in glob(path.format('*', '*')):
            ID = fname.split('/')[-1].split('.')[0]
            distribution, interaction  = ID.split('_')
            interactions.add(interaction)
            if distribution not in stats:
                stats[distribution] = {}
            stats[distribution][interaction] = DensityEstimate.read(fname)

        self = cls([], list(interactions), {}, {})
        self.stats = stats
        return self

    @classmethod
    def read_proteins(cls, path):
        stats = {}
        for fname in glob(path.format('*', '*', '*')):
            ID = fname.split('/')[-1].split('.')[0]
            protein, interaction, distribution = ID.split('-')
            if distribution not in stats:
                stats[distribution] = {}
            if interaction not in stats[distribution]:
                stats[distribution][interaction] = []
            stats[distribution][interaction] += [DensityEstimate.read(fname)]
        return stats

    def write_merged(self, merged):
        os.makedirs(merged, exist_ok=True)
        for dist, interactions in self.stats.items():
            for interaction, de in interactions.items():
                fname = '{}/{}_{}.txt'.format(merged, dist, interaction)
                with open(fname, 'w') as fp:
                    fp.write(str(de)+'\n')

    def plot_merged(self, plot):
        import matplotlib.pyplot as plt
        f, ax = plt.subplots(2, len(self.interactions),
                             figsize=(3*len(self.interactions), 3.5))

        for i, interaction in enumerate(self.interactions):
            nat = self.stats['native'][interaction]
            ref = self.stats['reference'][interaction]
            ax[0, i].set_title(interaction)
            ax[0, i].plot(nat.x, nat.fx, c='g')
            ax[0, i].plot(ref.x, ref.fx, c='k')
            ax[1, i].plot(nat.x, -np.log(nat.fx) + np.log(ref.fx), c='b')
            ax[0, i].set_ylim(0)
            ax[0, i].set_xlim(nat.x[0], nat.x[-1])
            ax[1, i].set_xlim(nat.x[0], nat.x[-1])
        ax[0, 0].set_ylabel('Frequency')
        ax[1, 0].set_ylabel('Energy')
        plt.savefig(plot)

    ############################################################################

    def _load(self):
        stats = {d: {i:[] for i in self.interactions}
                 for d in self.distributions}
        for protein in self.proteins:
            protein_stats = self._load_protein(protein)
            for d in self.distributions:
                for i in self.interactions:
                    stats[d][i] += [protein_stats[d][i]]
        return self._merge(stats)

    def _load_protein(self, protein):
        if self.path is not None:
            stats = self._read_protein(protein)
            if stats:
                return stats

        print('Computing stats for {}'.format(protein))
        
        prot = Protein(protein, self.settings, self.paths)
        ligands = prot.lm.get_xdocked_ligands(self.settings['n_ligs'])
        print(ligands)
        prot.load_docking(ligands, load_fp=True,
                          load_mcss='mcss' in self.interactions,
                          load_shape='shape' in self.interactions)

        stats = {d: {i: [] for i in self.interactions}
                 for d in self.distributions}
        for j, ligand1 in enumerate(ligands):
            for ligand2 in ligands[j+1:]:
                ligand_stats = self._load_ligand_pair(prot, protein, ligand1, ligand2)
                for d in self.distributions:
                    for i in self.interactions:
                        stats[d][i] += [ligand_stats[d][i]]
        stats = self._merge(stats)

        if self.path is not None:
            self._write_protein(protein, stats)
        return stats

    def _read_protein(self, protein):
        stats = {}
        for distribution in self.distributions:
            stats[distribution] = {}
            for interaction in self.interactions:
                fname = self.path.format(protein, interaction, distribution)
                try:
                    stats[distribution][interaction] = DensityEstimate.read(fname)
                except:
                    return {}
        return stats

    def _write_protein(self, protein, stats):
        if self.path is None:
            return

        for distribution, interactions in stats.items():
            for interaction, density_estimate in interactions.items():
                fname = self.path.format(protein, interaction, distribution)
                density_estimate.write(fname)

    def _load_ligand_pair(self, prot, protein, ligand1, ligand2):
        lig_pair = LigPair(prot.docking[ligand1],
                           prot.docking[ligand2],
                           self.interactions,
                           prot.lm.mcss  if 'mcss' in self.interactions else None,
                           prot.lm.shape if 'shape' in self.interactions else None,
                           self.settings['max_poses'])

        stats = {d: {} for d in self.distributions}
        for interaction in self.interactions:
            X_native, X_ref = self._get_interaction_scores(lig_pair, interaction)
            if 'mcss_domain' not in self.settings:
                self.settings['mcss_domain'] = (0, 15)
            domain = self.settings['mcss_domain'] if interaction == 'mcss' else (0, 1)
            sd = self.settings['stats_sd']*(domain[1]-domain[0])

            stats['native'][interaction] = DensityEstimate(domain=domain, sd=sd, reflect=True)
            stats['reference'][interaction] = DensityEstimate(domain=domain, sd=sd, reflect=True)
            stats['native'][interaction].fit(X_native)
            stats['reference'][interaction].fit(X_ref)
        return stats

    def _get_interaction_scores(self, lig_pair, interaction):
        n1 = min(len(lig_pair.l1.poses), self.settings['max_poses'])
        n2 = min(len(lig_pair.l2.poses), self.settings['max_poses'])
        
        X_native = []
        for r1 in range(n1):
            for r2 in range(n2):
                pp_x = lig_pair.get_feature(interaction, r1, r2)
                if pp_x is not None and lig_pair.correct(r1, r2):
                    X_native += [pp_x]
        X_ref = []
        for r1 in range(n1):
            for r2 in range(n2):
                pp_x = lig_pair.get_feature(interaction, r1, r2)
                if pp_x is not None:
                    X_ref += [pp_x]
        return np.array(X_native), np.array(X_ref)

    def _merge(self, stats):
        merged = {}
        for d, interactions in stats.items():
            merged[d] = {}
            for i, des in interactions.items():
                if des:
                    merged[d][i] = DensityEstimate.merge(des)
        return merged

def compute(params, paths, feature_defs, path, proteins, merged=None):
    path += '/{}-{}-{}.de'
    stats = Statistics(proteins, feature_defs, params, paths, path=path)

    if merged is not None:
        stats.write_merged(merged)

def plot(merged):
    stats = Statistics.read_merged(merged + '/{}_{}.txt')
    stats.plot_merged('{}/stats.pdf'.format(merged))

import numpy as np
from score.density_estimate import DensityEstimate
from score.pairs import LigPair
from glob import glob

class Statistics:
    def __init__(self, proteins, interactions, settings, paths, path=None):
        self.proteins = proteins
        self.interactions = interactions
        self.settings = settings
        self.distributions = ['native', 'reference']
        self.paths = paths
        self.path = path

        if 'reference_poses' not in settings:
            settings['reference_poses'] = settings['max_poses']
        if 'native_poses' not in settings:
            settings['native_poses'] = settings['max_poses']

        assert settings['max_poses'] >= settings['reference_poses']
        assert settings['max_poses'] >= settings['native_poses']

        self.stats = self._load()

    @classmethod
    def read(cls, path, paths, ligands_equal, considered_proteins=None):
        proteins = set()
        interactions = set()
        settings = {'ligands_equal': ligands_equal, 'max_poses': 100}
        stats = {}
        for fname in glob(path.format('*', '*', '*')):
            ID = fname.split('/')[-1].split('.')[0]
            protein, interaction, distribution = ID.split('-')
            if considered_proteins and protein not in considered_proteins:
                continue
            proteins.add(protein)
            interactions.add(interaction)
            if distribution not in stats:
                stats[distribution] = {}
            if interaction not in stats[distribution]:
                stats[distribution][interaction] = []
            stats[distribution][interaction] += [DensityEstimate.read(fname)]

        self = cls(list(proteins), list(interactions), settings, paths, path=path)
        self.stats = self._merge(stats, ligands_equal)
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

    def likelihood(self, interaction, distribution, x):
        return self.stats[distribution][interaction](x)

    def likelihood_ratio(self, interaction, x):
        n = self.likelihood(interaction, 'native', x)
        r = self.likelihood(interaction, 'reference', x)
        return n / r

    def domain(self, interaction):
        n = self.stats['native'][interaction].x[[0, -1]]
        r = self.stats['reference'][interaction].x[[0, -1]]
        assert np.all(n == r)
        return n

    def trace(self, interaction, points=100):
        low, high = self.domain(interaction)
        x = np.linspace(low, high, points)
        n = self.likelihood(interaction, 'native', x)
        r = self.likelihood(interaction, 'reference', x)
        return x, n, r

    def data_loglikelihood(self, protein):
        loglikelihood = {d: {i:0 for i in self.interactions}
                         for d in self.distributions}
        prot = Protein(protein)
        ligands = prot.lm.get_xdocked_ligands(self.settings['n_ligs'])
        prot.load_docking(ligands, load_fp=True, load_mcss=True)
        for j, ligand1 in enumerate(ligands):
            for ligand2 in ligands[j+1:]:
                lig_pair = LigPair(prot.docking[prot.lm.st].ligands[ligand1],
                                   prot.docking[prot.lm.st].ligands[ligand2],
                                   self.interactions, prot.lm.mcss,
                                   self.settings['max_poses'], self.settings['metric'])
                for i in self.interactions:
                    X_native, X_ref = self._get_interaction_scores(lig_pair, i)
                    loglikelihood['native'][i] += self.stats['native'][i].data_loglikelihood(X_native)
                    loglikelihood['reference'][i] += self.stats['reference'][i].data_loglikelihood(X_native)
        
        if self.settings['ligands_equal']:
            # Note that poses_equal is not implemented.
            n = len(ligands) * (len(ligands) - 1) / 2
            for distribution, interactions in loglikelihood.items():
                for interaction, _loglikelihood in interaction.items():
                    loglikelihood[distribution][interaction] = _loglikelihood / n
        return loglikelihood

    ############################################################################

    def _load(self):
        stats = {d: {i:[] for i in self.interactions}
                 for d in self.distributions}
        for protein in self.proteins:
            protein_stats = self._load_protein(protein)
            for d in self.distributions:
                for i in self.interactions:
                    stats[d][i] += [protein_stats[d][i]]
        return self._merge(stats, self.settings['ligands_equal'])

    def _load_protein(self, protein):
        if self.path is not None:
            stats = self._read_protein(protein)
            if stats:
                return stats

        print('Computing stats for {}'.format(protein))
        from containers import Protein
        prot = Protein(protein, self.settings, self.paths)
        ligands = prot.lm.get_xdocked_ligands(self.settings['n_ligs'])
        print(ligands)
        prot.load_docking(ligands, load_fp=True, load_mcss='mcss' in self.interactions)

        stats = {d: {i: [] for i in self.interactions}
                 for d in self.distributions}
        for j, ligand1 in enumerate(ligands):
            for ligand2 in ligands[j+1:]:
                ligand_stats = self._load_ligand_pair(prot, protein, ligand1, ligand2)
                for d in self.distributions:
                    for i in self.interactions:
                        stats[d][i] += [ligand_stats[d][i]]
        stats = self._merge(stats, self.settings['poses_equal'])

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
        lig_pair = LigPair(prot.docking[prot.lm.st][ligand1],
                           prot.docking[prot.lm.st][ligand2],
                           self.interactions,
                           prot.lm.mcss if 'mcss' in self.interactions else None,
                           self.settings['max_poses'],
                           self.settings['metric'])

        stats = {d: {} for d in self.distributions}
        for interaction in self.interactions:
            X_native, X_ref = self._get_interaction_scores(lig_pair, interaction)
            domain = (0, 15) if interaction == 'mcss' else (0, 1)
            sd = self.settings['stats_sd']*(domain[1]-domain[0])

            stats['native'][interaction] = DensityEstimate(domain=domain, sd=sd, reflect=True)
            stats['reference'][interaction] = DensityEstimate(domain=domain, sd=sd, reflect=True)
            stats['native'][interaction].fit(X_native)
            stats['reference'][interaction].fit(X_ref)
        return stats

    def _get_interaction_scores(self, lig_pair, interaction):
        X_native = []
        for (r1,r2), pp in lig_pair.pose_pairs.items():
            if max(r1, r2) >= self.settings['native_poses']: continue
            pp_x = lig_pair.get_feature(interaction, r1, r2)
            if pp_x is not None and pp.correct():
                X_native += [pp_x]
        X_ref = []
        for (r1,r2), pp in lig_pair.pose_pairs.items():
            if max(r1, r2) >= self.settings['reference_poses']: continue
            pp_x = lig_pair.get_feature(interaction, r1, r2)
            if pp_x is not None:
                X_ref += [pp_x]
        return np.array(X_native), np.array(X_ref)

    def _merge(self, stats, weight):
        merged = {}
        for d, interactions in stats.items():
            merged[d] = {}
            for i, des in interactions.items():
                merged[d][i] = DensityEstimate.merge(des, weight)
        return merged

def main(args):
    from settings import feature_defs, paths, stats
    import sys
    version, protein, path = args
    path += '/{}-{}-{}.de'
    Statistics([protein], feature_defs, stats[version], paths, path=path)


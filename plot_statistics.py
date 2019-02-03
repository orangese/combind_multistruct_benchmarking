import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import sys
# sys.path.append('../../dock/')
# sys.path.append('../../ifp/')
# sys.path.append('../../mcss/')
# sys.path.append('../../score/')
from score.density_estimate import DensityEstimate
from shared_paths import proteins, shared_paths, feature_defs
from score.statistics import statistics
import os
from glob import glob
# matplotlib inline

names = {
    'hbond'  : 'Hydrogen Bond',
    'sb'     : 'Salt Bridge',
    'pipi'   : 'Aromatic',
    'contact': 'Contact',
    'mcss'   : 'MCSS',
}

def main():
    for interaction in feature_defs:
        if interaction == 'mcss': continue
        # fig, ax = plt.subplots(figsize = (2.5, 2.0), dpi = 300, sharex=True, sharey=True)
        fig, ax = plt.subplots()
        m = 0
        # for n, protein in enumerate(proteins):
        for n, protein in enumerate(["B1AR"]):
            if protein in ['SIGMAR1', 'ERA', 'PYGM', 'CHRM3', 'TRPV1']: continue

            native_fname_format = '{}/{}/stats/{}/*-*-{}-native.de'.format(shared_paths['read_data'], protein,
                                                                          shared_paths['stats']['version'],
                                                                          interaction)
            reference_fname_format = '{}/{}/stats/{}/*-*-{}-reference.de'.format(shared_paths['read_data'], protein,
                                                                          shared_paths['stats']['version'],
                                                                          interaction)
            dude_fname_format = '{}/{}/stats/{}/*-*-{}-dude.de'.format(shared_paths['write_data'], protein,
                                                                          shared_paths['stats']['version'],
                                                                          interaction)
            
            # for fname in glob(native_fname_format):
            #     native = DensityEstimate.read(fname)
            #     ax.plot(native.x, native.fx, c = 'b', alpha = 0.05, lw=1)
            # for fname in glob(reference_fname_format):
            #     reference = DensityEstimate.read(fname)
            #     ax.plot(reference.x, reference.fx, c = 'g', alpha = 0.05, lw=1)
            # for fname in glob(dude_fname_format):
            #     dude = DensityEstimate.read(fname)
            #     ax.plot(dude.x, dude.fx, c = 'r', alpha = 0.05, lw=1)

            native_fname_format = '{}/{}/stats/{}/{}-{}-native.de'.format(shared_paths['read_data'], protein,
                                                                shared_paths['stats']['version'],
                                                                protein, interaction)
            reference_fname_format = '{}/{}/stats/{}/{}-{}-reference.de'.format(shared_paths['read_data'], protein,
                                                                shared_paths['stats']['version'],
                                                                protein, interaction)
            dude_fname_format = '{}/{}/stats/{}/{}-{}-dude_native.de'.format(shared_paths['write_data'], protein,
                                                                shared_paths['stats']['version'],
                                                                protein, interaction)
            native = DensityEstimate.read(native_fname_format)
            reference = DensityEstimate.read(reference_fname_format)
            dude = DensityEstimate.read(dude_fname_format)
            ax.set_title(protein, va = 'top')
            ax.plot(native.x, native.fx, c = 'b', lw = 3, label='Native')
            ax.plot(reference.x, reference.fx, c = 'g', lw = 3, label='Reference')
            ax.plot(dude.x, dude.fx, c = 'r', lw = 3, label='DUDE-Native')
            
            m = max(m, native.fx.max(), dude.fx.max(), reference.fx.max())
        
        ax.spines['left'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.set_yticks([])
        ax.set_yticklabels([])
        ax.set_ylim(0, m)
        ax.legend()

        if interaction == 'mcss':
            pass
        else:
            ax.set_xticks([0, 1])

        ax.set_xticklabels([0, 1], fontsize = 8)
        plt.suptitle(interaction)
        plt.savefig("plots/alldistr-dude-native-{}.png".format(interaction))

if __name__ == '__main__':
    main()

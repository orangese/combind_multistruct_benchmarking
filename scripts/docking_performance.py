import pickle
import sys
import os
sys.path.append(os.environ['COMBINDHOME'])
from containers import Protein
import utils
import config
import pandas as pd
import numpy as np
import click

family = pd.api.types.CategoricalDtype(['GPCR', 'Ion Channel', 'Nuclear Receptor', 'Transporter',
                                        'Reductase', 'Peptidase', 'Other'],
                                        ordered = True)

families = {
    'GPCR': ['P41595','P29274','P07700','P07550','P21554','P48039','P08172',
             'Q9UBS5','Q14416','Q14832','P41594','Q99835'],
    'Ion Channel': ['P19491', 'P22756', 'Q05586-Q12879', 'P42264'],
    'Transporter': ['Q7K4Y6','P31645','P11166'],
    'Nuclear Receptor': ['P10275','Q96RI1','P03372','P04150','P08235','P11473'],
    'Peptidase': ['P56817','P00760','P00742','P00734','P00749'],
    'Reductase': ['P15121','P28845','P00374','P16083'],
    'Other': ['P58154','Q99720','P36897','P34913','P09960','P24941','P28482',
              'P08246','P00489','P18031','O60885','P07900','Q9Y233']
}

drugs = {'GPCR': 0.33,
         'Ion Channel': 0.18,
         'Nuclear Receptor': 0.16,
         'Other': 0.13+0.03,
         'Peptidase': 0.03,
         'Reductase': 0.07,
         'Transporter': 0.07}

def drug_average(family):
    assert 'protein' not in family.index.names
    assert 'ligand' not in family.index.names
    weights = np.array([drugs[family] for family in family.index.get_level_values('family')])
    weighted = family * weights.reshape(-1, 1)
    return weighted.sum()

@click.command()
@click.option('--mcss-thresh', default=1.0)
@click.option('--best', is_flag=True)
@click.argument('version')
def main(version, mcss_thresh, best):

    paths = {'CODE': os.path.dirname(os.path.realpath(__file__)),
             'DATA': '/oak/stanford/groups/rondror/projects/ligand-docking/combind_bpp/combind2020',
             'PDB': '{ROOT}/structures/pdb.csv'}
    paths.update(config.PATHS)
    paths = utils.resolve(paths)
    params = config.STATS[version]

    proteins = utils.get_proteins(paths, [])

    data = []
    for prot in proteins:
        protein = Protein(prot, params, paths)
        ligands = protein.lm.get_pdb()
        protein.load_docking(ligands)
        for name, ligand in protein.docking.items():
            rmsds = [pose.rmsd for pose in ligand.poses]
            if None in rmsds:
                print(prot, name)
                continue
            if rmsds:
                best_rmsd = min(rmsds)
                top_rmsd = rmsds[0]
            else:
                best_rmsd = float('inf')
                top_rmsd = float('inf')

            data += [[prot, name, top_rmsd, best_rmsd]]

    data = pd.DataFrame(data, columns=['protein', 'lig', 'glide_rmsd', 'best_rmsd'])

    with open('stats_data/mcss_sizes.pkl', 'rb') as fp:
        mcss = pd.DataFrame(pickle.load(fp).items(), columns=['lig', 'mcss']).set_index('lig')
        data = data.join(mcss, on='lig')
    data = data.loc[data.mcss <= mcss_thresh]

    reverse = {v:k for k, vs in families.items() for v in vs}
    data['family'] = data['protein'].apply(lambda x: reverse[x]).astype(family)
    data = data.set_index(['family', 'protein', 'lig']).sort_index()

    data = data[['glide_rmsd', 'best_rmsd']].astype(float)
    data['total'] = 1

    if best:
        data = data.loc[data.best_rmsd <= 2.05]
    
    print((data <= 2.05).sum(axis=0))
    print()

    print((data <= 2.05).mean(axis=0))
    print()

    print(drug_average((data <= 2.05).groupby(level=0).mean()))
    print()

    print((data <= 2.05).groupby(level=0).mean())


main()


        # if top <= 2.0:
        #   print(prot, name, top, best)
import pandas as pd
import numpy as np
import click
import sys


family = pd.api.types.CategoricalDtype(['GPCR', 'Ion Channel', 'Nuclear Receptor', 'Transporter',
                                        'Peptidase', 'Other'],
                                        ordered = True)

families = {
    'GPCR': ['5HT2B', 'A2AR', 'B1AR', 'B2AR', 'SMO', 'MGLUR5'],
    'Ion Channel': ['P19491', 'P22756', 'Q05586-Q12879'],
    'Transporter': ['SLC6A4', 'GLUT1', 'DAT'],
    'Nuclear Receptor': ['NR3C2', 'NR3C1', 'AR', 'VDR', 'ERA'],
    'Peptidase': ['F2', 'F10', 'F11', 'PLAU', 'P00760', 'BACE1'],
    'Other': ['CDK2', 'PYGM', 'PTPN1', 'BRD4', 'HSP90AA1', 'PDE10A', 'SIGMAR1', 'ELANE', 'DHFR']
}

drugs = {'GPCR': 0.33,
         'Ion Channel': 0.18,
         'Nuclear Receptor': 0.16,
         'Other': 0.20+0.03,
         'Peptidase': 0.03,
         'Transporter': 0.07}


def drug_average(family):
    assert 'protein' not in family.index.names
    assert 'ligand' not in family.index.names
    weights = np.array([drugs[family] for family in family.index.get_level_values('family')])
    weighted = family * weights.reshape(-1, 1)
    return weighted.sum()

@click.command()
@click.option('--details', is_flag=True)
@click.option('--xtal', is_flag=True)
@click.option('--protein', default=-7)
@click.argument('fnames', nargs=-1)
def main(fnames, details, xtal, protein):
    data = []
    for path in fnames:
        prot = path.split('/')[protein]
        _data = pd.read_csv(path)
        _data = _data.iloc[:-1]
        if xtal:
            _data = _data[1:]
        if not len(_data):
            print(path)

        _data['protein'] = prot
        data += [_data]
    data = pd.concat(data)

    reverse = {v:k for k, vs in families.items() for v in vs}
    data['family'] = data['protein'].apply(lambda x: reverse[x]).astype(family)
    data = data.set_index(['family', 'protein', 'lig']).sort_index()

    data = data[['combind_rmsd', 'glide_rmsd', 'best_rmsd']].astype(float)
    data['total'] = 1

    print((data <= 2.05).sum(axis=0))
    print()

    print((data <= 2.05).groupby(level=0).mean())
    print()
    print(drug_average((data <= 2.05).groupby(level=0).mean()))

    if details:
        combind = data['combind_rmsd'] <= 2.05
        glide = data['glide_rmsd'] <= 2.05
        best = data['best_rmsd'] <= 2.05
        
        print('ComBind and not Glide')
        print(data[combind&~glide])
        print()

        print('Glide and not ComBind')
        print(data[~combind&glide])
        print()

        print('Best and not ComBind')
        print(data[~combind&best])
        print()


main()

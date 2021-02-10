import pandas as pd
import numpy as np
import click
import sys
import pickle

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
@click.option('--mcss-thresh', default=1.0)
@click.option('--best', is_flag=True)
@click.argument('fnames', nargs=-1)
def main(fnames, details, xtal, protein, mcss_thresh, best):
    data = []
    for path in fnames:
        prot = path.split('/')[protein]
        _data = pd.read_csv(path)
        _data = _data.iloc[:-1]
        _data = _data.loc[_data.COMBIND_RMSD!='None']
        _data['ID'] = [x.split('-to-')[0] for x in _data['ID']]
        if xtal:
            _data = _data[1:]
        if not len(_data):
            print(path)

        _data['protein'] = prot
        data += [_data]
    data = pd.concat(data)

    with open('/home/users/jpaggi/bpp_data_mcss.pkl', 'rb') as fp:
        mcss = pd.DataFrame(pickle.load(fp).items(), columns=['ID', 'mcss']).set_index('ID')

    with open('/home/users/jpaggi/ionchannels_mcss.pkl', 'rb') as fp:
        mcss = pd.concat([mcss, pd.DataFrame(pickle.load(fp).items(), columns=['ID', 'mcss']).set_index('ID')])
    
    data = data.join(mcss, on='ID')
    data = data.loc[data.mcss < mcss_thresh]

    reverse = {v:k for k, vs in families.items() for v in vs}
    data['family'] = data['protein'].apply(lambda x: reverse[x]).astype(family)
    data = data.set_index(['family', 'protein', 'ID']).sort_index()

    data = data[['COMBIND_RMSD', 'GLIDE_RMSD', 'BEST_RMSD']].astype(float)
    data['total'] = 1

    if best:
        data = data.loc[data.BEST_RMSD < 2.05]

    print((data <= 2.05).sum(axis=0))
    print()

    print((data <= 2.05).mean(axis=0))
    print()

    print(drug_average((data <= 2.05).groupby(level=0).mean()))
    print()

    print((data <= 2.05).groupby(level=0).mean())

    if details:
        combind = data['COMBIND_RMSD'] <= 2.05
        glide = data['GLIDE_RMSD'] <= 2.05
        best = data['BEST_RMSD'] <= 2.05
        
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

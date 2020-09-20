import pandas as pd
import numpy as np

def read_ifp(csv):
    df = pd.read_csv(csv)
    df.loc[df.label=='hbond_acceptor', 'protein_res'] = \
        [res+'acceptor' for res in df.loc[df.label=='hbond_acceptor', 'protein_res']]
    df.loc[df.label=='hbond_donor', 'protein_res'] = \
        [res+'donor' for res in df.loc[df.label=='hbond_donor', 'protein_res']]
    df.loc[df.label=='hbond_acceptor', 'label'] = 'hbond'
    df.loc[df.label=='hbond_donor', 'label'] = 'hbond'

    df = df.set_index(['label', 'pose',  'protein_res'])
    df = df.sort_index()
    return df

def ifp_tanimoto(ifp1, ifp2, feature):
    ifp1 = read_ifp(ifp1)
    ifp2 = read_ifp(ifp2)
    n1 = max(ifp1.index.get_level_values(level=1))+1
    n2 = max(ifp2.index.get_level_values(level=1))+1

    tanimotos = np.zeros((n1, n2))+0.5
    if feature not in ifp1.index.get_level_values(0):
        return tanimotos
    for i, _ifp1 in ifp1.loc[feature].groupby(['pose']):
        _ifp1 = _ifp1.droplevel(0)
        if feature not in ifp2.index.get_level_values(0):
            continue
        for j, _ifp2 in ifp2.loc[feature].groupby(['pose']):
            _ifp2 = _ifp2.droplevel(0)

            joined = _ifp1.join(_ifp2, how='outer', lsuffix='1', rsuffix='2')
            joined = joined.fillna(0)

            overlap = sum((joined['score1']*joined['score2'])**0.5)
            total = sum(joined['score1']+joined['score2'])

            tanimotos[i, j] = (1 + overlap) / (2 + total - overlap)
    return tanimotos

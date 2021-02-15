import pandas as pd
import numpy as np

def read_ifp(csv):
    """
    Reads IFP file and merges hbond acceptors and donors.

    Setting the label to hbond for the hbond_donors and hbond_acceptors while
    changing the residue names allows for only donor+donor or acceptor+acceptor
    to be counted as overlapping, but them to be merged into the same similarity
    measure.
    """
    df = pd.read_csv(csv)

    mask = df.label=='hbond_acceptor'
    df.loc[mask, 'protein_res'] = [res+'acceptor' for res in df.loc[mask, 'protein_res']]
    df.loc[mask, 'label'] = 'hbond'
    
    mask = df.label=='hbond_donor'
    df.loc[mask, 'protein_res'] = [res+'donor' for res in df.loc[mask, 'protein_res']]
    df.loc[mask, 'label'] = 'hbond'
    return df

def ifp_tanimoto(ifp1, ifp2, feature):
    """
    Computes the tanimoto distance between ifp1 and ifp2 for feature.
    """
    ifp1 = read_ifp(ifp1)
    ifp2 = read_ifp(ifp2)
    n1 = max(ifp1.pose)+1
    n2 = max(ifp2.pose)+1

    ifp1 = ifp1.loc[ifp1.label == feature]
    ifp2 = ifp2.loc[ifp2.label == feature]

    interactions = sorted(set(ifp1.protein_res).union(set(ifp2.protein_res)))

    X1 = np.zeros((n1, 1, len(interactions)))
    for i, row in ifp1.iterrows():
        X1[row.pose, 0, interactions.index(row.protein_res)] = row.score

    X2 = np.zeros((1, n2, len(interactions)))
    for i, row in ifp2.iterrows():
        X2[0, row.pose, interactions.index(row.protein_res)] = row.score

    overlap = np.sqrt(X1*X2).sum(axis=2)
    total = X1.sum(axis=2) + X2.sum(axis=2)

    return (1 + overlap) / (2 + total - overlap)

import os
import sys

def load_features(fpath=None):
    txt_file = 'ligands/ligand_features.txt'
    if fpath is not None:
        txt_file = fpath + '/' + txt_file
        
    if not os.path.exists(txt_file): 
        return {}

    lig_features = {}
    with open(txt_file) as f:
        for line in f:
            line = line.strip().split(':')
            if len(line) != 2: continue
            lig_features[line[0]] = line[1].split(',')

    return lig_features

def output_features():
    from schrodinger.structure import StructureReader

    done = load_features()
    not_done = {}
    for lig in os.listdir('ligands/unique'):

        if lig[-3:] != 'mae': continue
        if lig[0] == '.': continue

        lig = lig.split('.')[0]
        if lig in done: continue

        not_done[lig] = []
        st = StructureReader('ligands/unique/{}.mae'.format(lig)).next()

        for a in st.atom:
            if a.formal_charge != 0:
                not_done[lig].append('chrg')
                break

        if len([ri for ri in st.ring if ri.isAromatic() or ri.isHeteroaromatic()]) > 0:
            not_done[lig].append('ring')

    if len(not_done) > 0:
        with open('ligands/ligand_features.txt', 'a') as f:
            for lig, feats in sorted(not_done.items()):
                f.write('{}:{}\n'.format(lig,','.join(feats)))




import click
import pandas as pd
import os
from schrodinger.structure import SmilesStructure
from schrodinger.structutils.analyze import find_common_substructure

def get_structure(smiles):
    return SmilesStructure(smiles).get2dStructure()

def heavy_atoms(st):
    return sum(atom.element != 'H' for atom in st.atom)

def mcss_size(st1, st2, size_method=min):
    mcss = find_common_substructure([st1, st2], atomTyping=12)
    mcss_size = len(mcss[0][0])
    ligand_size = size_method(heavy_atoms(st1), heavy_atoms(st2))
    return mcss_size / ligand_size

def unique(ligands_st, query_st):
    for ligand_st in ligands_st:
        if ligand_st.isEquivalent(query_st):
            return False
    return True

def diverse(ligands_st, query_st):
    for ligand_st in ligands_st:
        if mcss_size(ligand_st, query_st) > 0.8:
            return False
    return True

def pick_affinity(query, ligands, n_helpers):
    ligands = sorted(ligands, key=lambda x: x[2])
    query_st = get_structure(query[1])
    helpers, helpers_st = [], []
    for ligand in ligands:
        ligand_st = get_structure(ligand[1])
        if unique(helpers_st+[query_st], ligand_st):
            helpers += [ligand]
            helpers_st += [ligand_st]
        if len(helpers) == n_helpers:
            break
    return helpers

def pick_affinity_diverse(query, ligands, n_helpers):
    ligands = sorted(ligands, key=lambda x: x[2])
    query_st = get_structure(query[1])
    helpers, helpers_st = [], []
    for ligand in ligands:
        ligand_st = get_structure(ligand[1])
        if unique(helpers_st+[query_st], ligand_st) and diverse(helpers_st, ligand_st):
            helpers += [ligand]
            helpers_st += [ligand_st]
        if len(helpers) == n_helpers:
            break
        print(len(helpers))
    return helpers

def pick_mcss(query, ligands, n_helpers):
    query_st = get_structure(query[1])
    helpers, helpers_st = [], []
    for ligand in ligands:
        ligand_st = get_structure(ligand[1])
        if unique(helpers_st+[query_st], ligand_st):
            helpers += [ligand]
            helpers_st += [ligand_st]

    mcss = {}
    for ligand, st in zip(helpers, helpers_st):
        mcss[ligand[0]] = mcss_size(query_st, st)

    helpers = sorted(helpers, key=lambda x: x[2])
    helpers = sorted(helpers, key=lambda x: mcss[x[0]])
    return helpers[:n_helpers]

@click.command()
@click.argument('queries_csv')
@click.argument('helpers_csv')
@click.argument('out_csv_dir')
@click.option('--criteria', default='affinity')
@click.option('--n_helpers', default=25)
@click.option('--mcss-version', default='mcss16')
@click.option('--mcss-size', default='min', type=click.Choice(['min', 'max', 'avg']))
@click.option('--mcss-diverse-cut', default=0.8)
@click.option('--affinity-cut', default=1000)
def main(queries_csv, helpers_csv, out_csv_dir, criteria, n_helpers,
         mcss_version, mcss_size, mcss_diverse_cut, affinity_cut):
    queries = pd.read_csv(queries_csv)
    helpers = pd.read_csv(helpers_csv)
    helpers = helpers.loc[helpers['AFFINITY'] < affinity_cut]

    queries = [v.to_list() for _, v in queries[['ID', 'SMILES']].iterrows()]
    helpers = [v.to_list() for _, v in helpers[['ID', 'SMILES', 'AFFINITY']].iterrows()]

    for query in queries:
        fname = '{}/{}-{}.csv'.format(out_csv_dir, query[0], criteria)
        if os.path.exists(fname):
            continue
        print(fname)
        
        if criteria == 'affinity':
            _helpers = pick_affinity(query, helpers, n_helpers)
        elif criteria == 'affinity_diverse':
            _helpers = pick_affinity_diverse(query, helpers, n_helpers)
        elif criteria == 'mcss':
            _helpers = pick_mcss(query, helpers, n_helpers)
        else:
            assert False

        with open(fname, 'w') as fp:
            fp.write('ID,SMILES,AFFINITY\n')
            fp.write(','.join(map(str, query + [0])) + '\n')
            for helper in _helpers:
                fp.write(','.join(map(str, helper)) + '\n')

main()

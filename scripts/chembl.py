import sqlite3
import os
import pandas as pd
from schrodinger.structure import SmilesStructure, StructureReader
from schrodinger.structutils.analyze import generate_smiles
import click
import tempfile

MACROCYCLE_THRESH = 8
AFFINITY_THRESH = 1000
MOLW_THRESH = 800

class CHEMBLDB:
    def __init__(self, chembldb, uniprot_chembl):
        self.conn = sqlite3.connect('file:{}?mode=ro'.format(chembldb), uri=True)
        self.cur = self.conn.cursor()
        self.uniprot_chembl = pd.read_csv(uniprot_chembl, sep='\t', index_col=0)
        self.uniprot_chembl = self.uniprot_chembl['Cross-reference (ChEMBL)'].apply(lambda x: x.strip(';').split(';'))

    def __enter__ (self):
        return self
    def __exit__ (self, *exc):
        self.conn.close()

    def chembl_to_tid(self, chembl):
        self.cur.execute("SELECT tid FROM target_dictionary WHERE chembl_id=?", (chembl,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def tid_to_target_type(self, tid):
        self.cur.execute("SELECT target_type FROM target_dictionary WHERE tid=?", (tid,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def tid_to_assays(self, tid, protein_complex):
        if protein_complex:
            self.cur.execute("SELECT assay_id FROM assays WHERE tid=? AND (confidence_score=7 OR confidence_score=9)", (tid,))
        else:
            self.cur.execute("SELECT assay_id FROM assays WHERE tid=? AND confidence_score=9", (tid,))
        return [row[0] for row in self.cur.fetchall()]

    def assay_to_molregnos(self, assay):
        self.cur.execute("SELECT molregno FROM activities WHERE assay_id=?", (assay,))
        return [row[0] for row in self.cur.fetchall()]

    def molregno_to_smiles(self, molregno):
        self.cur.execute("SELECT canonical_smiles FROM compound_structures WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def molregno_to_molw(self, molregno):
        self.cur.execute("SELECT mw_freebase FROM compound_properties WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        if not len(rows):
            return 0
        return rows[0][0]

    def molregno_to_smiles(self, molregno):
        self.cur.execute("SELECT canonical_smiles FROM compound_structures WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        if not rows:
            return None
        return rows[0][0]

    def molregno_to_chemblid(self, molregno):
        self.cur.execute("SELECT chembl_id FROM molecule_dictionary WHERE molregno=?", (molregno,))
        rows = self.cur.fetchall()
        assert len(rows) == 1, rows
        return rows[0][0]

    def molregno_and_assay_to_activities(self, molregno, assay):
        self.cur.execute("SELECT standard_type, standard_value, standard_units, relation FROM activities WHERE molregno=? AND assay_id=?", (molregno, assay))
        return self.cur.fetchall()

    def chembl_to_activities(self, chembl, protein_complex):
        # chemblID, SMILES, MOLW, affinity
        activities = []
        tid = self.chembl_to_tid(chembl)
        for assay in self.tid_to_assays(tid, protein_complex):
            for molregno in self.assay_to_molregnos(assay):
                molw = self.molregno_to_molw(molregno)
                smiles = self.molregno_to_smiles(molregno)
                chembl_id = self.molregno_to_chemblid(molregno)
                for activity in self.molregno_and_assay_to_activities(molregno, assay):
                    activities += [[chembl_id, molw, smiles] + list(activity)]
        return pd.DataFrame(activities,
                            columns=['ligand_chembl_id', 'mw_freebase', 'canonical_smiles',
                                     'standard_type', 'standard_value',
                                     'standard_units', 'relation'])

    def uniprot_to_chembl(self, uniprot):
        for chembl_id in self.uniprot_chembl.loc[uniprot]:
            tid = self.chembl_to_tid(chembl_id)
            target_type = self.tid_to_target_type(tid)
            if target_type == 'SINGLE PROTEIN':
                return chembl_id

################################################################################

def get_chembl(uniprot, chembldb, uniprot_chembl):
    with CHEMBLDB(chembldb, uniprot_chembl) as chembldb:
        chembl = chembldb.uniprot_to_chembl(uniprot)
    return chembl

def get_activities(chembl, chembldb, uniprot_chembl, protein_complex):
    with CHEMBLDB(chembldb, uniprot_chembl) as chembldb:
        activities = chembldb.chembl_to_activities(chembl, protein_complex)
    activities['target_chembl_id'] = chembl
    return activities

def filter_activities(activities, activity_type):
    # Standardize units
    m = {'M': 10**9, 'mM': 10**6, 'uM': 10**3, 'pM': 10**-3}
    for unit, relation in m.items():
        mask = activities['standard_units'] == unit
        activities.loc[mask, 'standard_value'] *= relation
        activities.loc[mask, 'standard_units'] = unit

    # Filter
    mask = activities['standard_value'].notna()
    print('Removing {} rows b/c standard_value is na'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = activities['standard_value'] != 0
    print('Removing {} rows b/c standard_value is 0'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = activities['canonical_smiles'] != None
    print('Removing {} rows b/c canonical_smiles is None'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = activities['mw_freebase'] <= MOLW_THRESH+100 # Not desalted yet.
    print('Removing {} rows b/c mw_freebase > {}'.format(len(mask)-sum(mask),
                                                         MOLW_THRESH+100))
    activities = activities.loc[mask]

    mask = activities['standard_units'] == 'nM'
    print('Removing {} rows b/c standard_units != nM'.format(len(mask)-sum(mask)))
    print('Set of offending values is {}'.format(set(activities[~mask]['standard_units'])))
    activities = activities.loc[mask]

    mask = activities['relation'] == '='
    print('Removing {} rows b/c relation != ='.format(len(mask)-sum(mask)))
    print('Set of offending values is {}'.format(set(activities[~mask]['relation'])))
    activities = activities.loc[mask]

    if activity_type == 'all':
        activity_types = ['EC50', 'IC50', 'Ki', 'Kd']
    else:
        activity_types = [activity_type]
    mask = activities['standard_type'].isin(activity_types)
    print('Removing {} rows b/c standard_type not in {}'.format(len(mask)-sum(mask),
                                                                 activity_types))
    print('Set of offending values is {}'.format(set(activities[~mask]['standard_type'])))
    activities = activities.loc[mask]

    # Average
    keys = ['ligand_chembl_id', 'standard_type']
    averages = activities.loc[:, keys+['standard_value']].groupby(keys).mean()
    activities = activities.groupby(keys, as_index=False).first()
    activities['standard_value'] = [averages.loc[tuple([row[key] for key in keys])]['standard_value']
                                    for _, row in activities.iterrows()]
    return activities

################################################################################

def get_structure(smiles):
    smi = SmilesStructure(smiles)
    try:
        st = smi.get3dStructure(True)
        stereo = True
    except:
        try:
            st = smi.get3dStructure(False)
            stereo = False
        except:
            print('Error processing {}'.format(smiles))
            return SmilesStructure('C').get3dStructure(), False

    # Desalt.
    atoms = []
    for mol in st.molecule:
        if len(mol.atom) > len(atoms):
            atoms = [a.index for a in mol.atom]
    st = st.extract(atoms)
    return st, stereo

def is_macrocycle(st):
    ring_sizes = [0]+[len(ring.atom) for ring in st.ring]
    return max(ring_sizes) > MACROCYCLE_THRESH

def _get_properties(smiles):
    properties = {}
    st, properties['stereo'] = get_structure(smiles)
    properties['SMILES'] = generate_smiles(st)
    properties['macrocycle'] = is_macrocycle(st)
    properties['molw'] = st.total_weight
    return pd.Series(properties)

def get_properties(activities):
    properties = activities.canonical_smiles.apply(_get_properties)
    return pd.concat([activities, properties], axis=1)

def filter_properties(activities):
    mask = activities['molw'] <= MOLW_THRESH
    print('Removing {} rows b/c molw > {}'.format(len(mask)-sum(mask),
                                                  MOLW_THRESH))
    activities = activities.loc[mask]

    mask = ~activities['macrocycle']
    print('Removing {} rows b/c macrocycle'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]

    mask = activities['stereo']
    print('Removing {} rows b/c ambiguous stereochemistry'.format(len(mask)-sum(mask)))
    activities = activities.loc[mask]
    return activities


@click.command()
@click.option('--protein-complex', is_flag=True)
@click.option('--activity_type', default='all')
@click.argument('uniprot_or_chembl')
@click.argument('chembldb', default='/oak/stanford/groups/rondror/users/jpaggi/pldb_data/raw/chembl_25.db')
@click.argument('uniprot_chembl', default='/oak/stanford/groups/rondror/users/jpaggi/pldb_data/raw/uniprot-chembl.tsv')
def main(protein_complex, activity_type, uniprot_or_chembl,
         chembldb, uniprot_chembl):
    if uniprot_or_chembl[:6] == 'CHEMBL':
        chembl = uniprot_or_chembl
    else:
        chembl = get_chembl(uniprot_or_chembl, chembldb, uniprot_chembl)
    if chembl is None:
        exit()
    activities = get_activities(chembl, chembldb, uniprot_chembl, protein_complex)
    print(chembl, len(activities))
    activities = filter_activities(activities, activity_type)
    activities = get_properties(activities)
    activities = filter_properties(activities)
    activities['AFFINITY'] = activities['standard_value']
    activities['ID'] = activities['ligand_chembl_id']
    activities.to_csv('{}_{}.csv'.format(chembl, activity_type), index=False)

if __name__ == '__main__':
    main()

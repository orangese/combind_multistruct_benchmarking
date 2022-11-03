#!/bin/env python

import pandas as pd
import numpy as np
import click
import os
from glob import glob
from schrodinger.structure import StructureReader, StructureWriter
from utils import *

###############################################################################

# Defaults
stats_root = os.environ['COMBINDHOME']+'/stats_data/default'
mcss_version = 'mcss16'
shape_version = 'pharm_max'
ifp_version = 'rd1'

@click.group()
def main():
    pass



import os
from schrodinger.structure import StructureReader

def split_complex(st, pdb_id, template):
    os.system('mkdir -p structures/proteins structures/ligands')
    lig_path = 'structures/ligands/{}_lig_to_{}.mae'.format(pdb_id, template)
    prot_path = 'structures/proteins/{}_prot_to_{}.mae'.format(pdb_id, template)

    if not os.path.exists(lig_path) and len([a.index for a in st.atom if a.chain == 'L']) > 0:
        lig_st = st.extract([a.index for a in st.atom if a.chain == 'L'])
        lig_st.title = '{}_lig_to_{}'.format(pdb_id, template)
        lig_st.write(lig_path)
    
    if not os.path.exists(prot_path):
        prot_st = st.extract([a.index for a in st.atom if a.chain != 'L'])
        prot_st.title = '{}_prot_to_{}'.format(pdb_id, template)
        prot_st.write(prot_path)

def struct_sort(structs):
    for struct in structs:
        for template in structs:
            opt_complex = 'structures/aligned/{}/rot-{}_query_to_{}.mae'.format(struct, struct, struct_template)

            if os.path.exists(opt_complex):
                comp_st = next(StructureReader(opt_complex))
                split_complex(comp_st, struct, template)




@main.command()
@click.argument('docking', default='docking/*/*_pv.maegz')
@click.argument('crystal', default='structures/ligands/*_lig_to_*.mae')
def rmsd_all(docking, crystal):
    """
    Compute rmsd of docked poses to a reference pose.

    docking and crystal should be paths to the docked poses in pose viewer
    format and the reference poses. Non-expanded patterns capturing multiple
    docking results and references (for different ligands) can be provided to
    process multiple files at a time.

    It is required that the docked poses and crystal ligands have the same
    name. For the docking the name is the part of the basename before '-to-' and
    for the reference poses it is the part of the basename before '_lig'.
    """
    from dock.dock import rmsd

    docking = glob(docking)
    crystal = glob(crystal)

    def docking_to_name(path):
        path = path.split('/')[-1]
        ligand = path.split('-to-')[0]
        ligand = ligand.split('_lig')[0]
        struct = path.split('-to-')[1]
        return ligand, struct

    def crystal_to_name(path):
        path = path.split('/')[-1]
        ligand = path.split('_lig_to_')[0]
        template = path.split('_lig_to_')[1]
        return ligand, template

    for docking_path in docking:
        # if os.path.exists(out):
            # continue
        
        name = docking_to_name(docking_path)
        crystal_path = [crystal_path for crystal_path in crystal
                        if crystal_to_name(crystal_path) == name]

        if len(crystal_path) == 0:
            print('No crystal pose for {}: {}'.format(name, docking_path))
            continue
        if len(crystal_path) > 1:
            print('Multiple crystal poses for {}: {}. Doing nothing.')
            continue

        crystal_path = crystal_path[0]
        print('Computing rmsd for {} to {}.'.format(docking_path, crystal_path))
        for path in crystal_path:
            out = docking_path.replace('.maegz', '_rmsd_to.npy')
            rmsds = rmsd(crystal_path, docking_path)
        np.save(out, rmsds) 
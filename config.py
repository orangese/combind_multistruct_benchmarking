from utils import StringFunction

STATS = {'default': {'version'         : 'default',
                     'ifp_version'     : 'default',
                     'mcss_version'    : 'mcss16',
                     'mcss_func'       : StringFunction('min'),
                     'mcss_rel_min'    : 0.5,
                     'mcss_abs_min'    : 10,
                     'mcss_domain'     : (0, 6),
                     'docking_version' : 'confgen_es4',
                     'pdb_order'       : 'First',
                     'native_thresh'   : 2.0,
                     'n_ligs'          : 20,
                     'max_poses'       : 100,
                     'stats_sd'        : 0.03,
                     'stats_points'    : 100}}

IFP = {'default': {'version'           : 'default',
                   'hbond_dist_opt'    : 2.5,
                   'hbond_dist_cut'    : 3.0,
                   'hbond_angle_opt'   : 60.0,
                   'hbond_angle_cut'   : 90.0,
                   'sb_dist_opt'       : 4.0,
                   'sb_dist_cut'       : 5.0,
                   'contact_scale_opt' : 1.25,
                   'contact_scale_cut' : 1.75},}

FEATURE_DEFS = {
    'mcss'           : [],
    'sb'             : [1],
    'hbond'          : [2,3],
    'hbond_donor'    : [2],
    'hbond_acceptor' : [3],
    'contact'        : [11]
}

PATHS = {'ROOT': '{DATA}/{protein}',

         'LIGANDS_ROOT': '{ROOT}/ligands',
         'LIGANDS': '{LIGANDS_ROOT}/{ligand}/{ligand}.mae',
         'PDB': '{ROOT}/structures/pdb.csv',

         'GRID_ROOT': '{ROOT}/docking/grids',
         'GRID': '{GRID_ROOT}/{struct}/{struct}.zip',
         'DOCK': '{ROOT}/docking/{docking_version}/{ligand}-to-{struct}',
         'DOCK_PV': '{DOCK}/{ligand}-to-{struct}_pv.maegz',
         'DOCK_REPT': '{DOCK}/{ligand}-to-{struct}.rept',
         'DOCK_RMSD': '{DOCK}/rmsd.csv',
         'IFP_ROOT':  '{ROOT}/ifp/{ifp_version}',
         'IFP':  '{IFP_ROOT}/{ligand}-to-{struct}-{docking_version}.fp',
         'MCSS': '{ROOT}/mcss/{mcss_version}',
         }

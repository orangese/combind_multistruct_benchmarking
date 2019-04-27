
stats = {'stats5': {'version'         : 'stats5',
                    'ifp_version'     : 'ifp3',
                    'mcss_version'    : 'mcss14',
                    'docking_version' : 'glide12',
                    'pdb_order'       : 'First',
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'weighting'       : 'absolute',
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,
                    'gscore_domain'   : (-16, 2),
                    'gscore_points'   : 1000,
                    'gscore_sd'       : 0.4},

         'stats6': {'version'         : 'stats6',
                    'ifp_version'     : 'ifp4',
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,
                    'weighting'       : 'unweighted',
                    'gscore_domain'   : (-16, 2),
                    'gscore_points'   : 1000,
                    'gscore_sd'       : 0.4},

         'stats7': {'version'         : 'stats7',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,
                    'weighting'       : 'unweighted',
                    'gscore_domain'   : (-16, 2),
                    'gscore_points'   : 1000,
                    'gscore_sd'       : 0.4},

         'stats8': {'version'         : 'stats8',
                    'ifp_version'     : 'ifp6',
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,
                    'weighting'       : 'unweighted',
                    'gscore_domain'   : (-16, 2),
                    'gscore_points'   : 1000,
                    'gscore_sd'       : 0.4},
         }

ifp = {'ifp3': {'version': 'ifp3',
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 3.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 5.0,
                'pipi_dist_opt': 4.5,
                'pipi_dist_cut': 6.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75},

       'ifp4': {'version'  : 'ifp4', # Cutoffs same as in ifp3, but ifp code changed.
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 3.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 5.0,
                'pipi_dist_opt': 4.5,
                'pipi_dist_cut': 6.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75},

       'ifp5': {'version'  : 'ifp5', # Cutoffs same as in ifp4, but ifp code changed.
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 3.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 5.0,
                'pipi_dist_opt': 4.5,
                'pipi_dist_cut': 6.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75},

       'ifp6': {'version'  : 'ifp6',
                'hbond_dist_opt': 2.75,
                'hbond_dist_cut': 5.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 3.0,
                'sb_dist_cut': 6.0,
                'pipi_dist_opt': 7.75,
                'pipi_dist_cut': 9.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75},
     }

feature_defs = {
    'mcss'           :[],
    'sb'             :[1],
    'hbond'          :[2,3],
    'hbond_donor'    :[2],
    'hbond_acceptor' :[3],
    'pipi'           :[6],
    'contact'        :[11]
}

shared_paths = { 
    'code'      : '~/combind',
    'data'      : "/scratch/PI/rondror/combind/bpp_data",
    'docking'   : 'confgen_es4',
    'mcss'      : 'mcss16',
    'stats'     : stats['stats6'],
    'ifp'       : ifp['ifp4'],
    'pdb_order' : 'First'
}

import os
exclude = []
proteins = [p for p in os.listdir(shared_paths['data']) if p[0] != '.' and p not in exclude]

assert shared_paths['ifp']['version'] == shared_paths['stats']['ifp_version']
assert shared_paths['docking'] == shared_paths['stats']['docking_version']
assert shared_paths['mcss'] == shared_paths['stats']['mcss_version']
assert shared_paths['pdb_order'] == shared_paths['stats']['pdb_order']

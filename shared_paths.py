
class StringFunction:
    """
    Wraps function such that it's string representation
    can be read with the eval() function.
    """
    def __init__(self, func):
        self.str = func
        self.func = eval(func)
    def __call__(self, *args):
        return self.func(*args)
    def __str__(self):
        return self.str
    def __repr__(self):
        return self.str

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

         'stats9': {'version'         : 'stats9',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : True,
                    'ligands_equal'   : True,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100},
         
         'stats10': {'version'        : 'stats10',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.03,
                    'stats_points'    : 100},

        'stats11': {'version'        : 'stats11',
                    'ifp_version'     : 'ifp7',
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100},

        'stats12': {'version'         : 'stats12', # Uses tanimoto coeff.
                    'ifp_version'     : 'ifp5',    # Must hardcode.
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
         'stats13': {'version'        : 'stats13', # Uses tanimoto coeff.
                    'ifp_version'     : 'ifp5',    # Must hardcode.
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : True,
                    'ligands_equal'   : True,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats14': {'version'         : 'stats14', # Uses tanimoto coeff.
                    'ifp_version'     : 'ifp5',    # Must hardcode.
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : True,
                    'ligands_equal'   : True,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats15': {'version'         : 'stats15', # Uses tanimoto coeff.
                    'ifp_version'     : 'ifp5',    # Must hardcode.
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : True,
                    'ligands_equal'   : True,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 1,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats16': {'version'         : 'stats16', # Uses psuedo-tanimoto coeff.
                    'ifp_version'     : 'ifp5',    # Must hardcode.
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : True,
                    'ligands_equal'   : True,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats17': {'version'         : 'stats17', # Uses psuedo-tanimoto coeff.
                    'ifp_version'     : 'ifp5',    # Must hardcode.
                    'mcss_version'    : 'mcss16',
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats18': {'version'         : 'stats18',
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
                    'stats_points'    : 100,},
        'stats19': {'version'         : 'stats19',
                    'metric'          : 'maxoverlap',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('max'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats20': {'version'         : 'stats20',
                    'metric'          : 'maxoverlap',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('max'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'XP',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats21': {'version'         : 'stats21',
                    'metric'          : 'tanimoto',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('min'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats22': {'version'         : 'stats22',
                    'metric'          : 'maxoverlap',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('min'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats23': {'version'         : 'stats23',
                    'metric'          : 'tanimoto',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('min'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.02,
                    'stats_points'    : 100,},
        'stats24': {'version'         : 'stats24', # min instead of geometric mean.
                    'metric'          : 'tanimoto',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('min'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
        'stats25': {'version'         : 'stats25', # Exclude interaction from stats if no poses make one.
                    'metric'          : 'tanimoto',
                    'ifp_version'     : 'ifp5',
                    'mcss_version'    : 'mcss16',
                    'mcss_func'       : StringFunction('min'),
                    'mcss_rel_min'    : 0.5,
                    'mcss_abs_min'    : 10,
                    'docking_version' : 'confgen_es4',
                    'pdb_order'       : 'First',
                    'poses_equal'     : False,
                    'ligands_equal'   : False,
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,},
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

       'ifp7': {'version'  : 'ifp7',
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 2.5,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 60.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 4.0,
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
    'data'      : "/oak/stanford/groups/rondror/users/jpaggi/bpp_data",
    'docking'   : 'confgen_es4',
    'mcss'      : 'mcss16',
    'stats'     : stats['stats21'],
    'ifp'       : ifp['ifp5'],
    'pdb_order' : 'First'
}

import os
exclude = ['D2']
proteins = [p for p in os.listdir(shared_paths['data']) if p[0] != '.' and p not in exclude]

assert shared_paths['ifp']['version'] == shared_paths['stats']['ifp_version']
assert shared_paths['docking'] == shared_paths['stats']['docking_version']
assert shared_paths['mcss'] == shared_paths['stats']['mcss_version']
assert shared_paths['pdb_order'] == shared_paths['stats']['pdb_order']

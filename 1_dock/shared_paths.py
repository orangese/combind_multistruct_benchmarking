stats = {'stats5': {'version'         : 'stats5',
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
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'weighting'       : 'unweighted',
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,
                    'gscore_domain'   : (-16, 2),
                    'gscore_points'   : 1000,
                    'gscore_sd'       : 0.4}
         }

ifp = {'ifp3': {'version'  : 'ifp3',
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 3.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 5.0,
                'pipi_dist_opt': 4.5,
                'pipi_dist_cut': 6.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75}
     }

# Cutoffs same as in ifp3, but ifp code changed.
ifp = {'ifp4': {'version'  : 'ifp4',
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 3.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 5.0,
                'pipi_dist_opt': 4.5,
                'pipi_dist_cut': 6.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75}
     }

# Increased cut distance so close interactions
# are written to file.
ifp = {'raw_ifp4': {'version'  : 'raw',
                'hbond_dist_opt': 2.5,
                'hbond_dist_cut': 5.0,
                'hbond_angle_opt': 60.0,
                'hbond_angle_cut': 90.0,
                'sb_dist_opt': 4.0,
                'sb_dist_cut': 6.0,
                'pipi_dist_opt': 4.5,
                'pipi_dist_cut': 9.0,
                'contact_scale_opt': 1.25,
                'contact_scale_cut': 1.75}
     }

feature_defs = {
    'mcss':[],
    'hbond':[2,3],
    'hbond_donor':[2],
    'hbond_acceptor':[3],
    'sb1':[0],
    'sb2':[1],
    'sb3':[4],
    'pipi':[6],
    'picat':[7,8],
    'contact':[11]
}

shared_paths = { 
    'code'      : '~/combind',
    'data'      : "/scratch/PI/rondror/combind/bpp_data/",
    'docking'   : 'glide12',
    'mcss'      : 'mcss14',
    'stats'     : stats['stats6'],
    'ifp'       : ifp['raw_ifp4']
}

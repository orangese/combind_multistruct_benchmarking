stats = {'stats5': {'version'         : 'stats5',
                    'native_thresh'   : 2.0,
                    'n_ligs'          : 20,
                    'max_poses'       : 100,
                    'weighting'       : 'absolute',
                    'stats_sd'        : 0.07,
                    'stats_points'    : 100,
                    'gscore_domain'   : (-16, 2),
                    'gscore_points'   : 1000,
                    'gscore_sd'       : 0.4}
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
    'data'      : "/Users/jpaggi/Documents/combind/combind_data/bpp_data/",
    'docking'   : 'glide12',
    'mcss'      : 'mcss14',
    'stats'     : stats['stats5'],
    'ifp'       : 'ifp3'
}

import os
import numpy as np
from glob import glob
from schrodinger.structure import StructureReader
from utils import *

IFP = {'rd1':    {'version'           : 'rd1',
                   'level'             : 'residue',
                   'hbond_dist_opt'    : 3.5,
                   'hbond_dist_cut'    : 4.0,
                   'hbond_angle_opt'   : 60.0,
                   'hbond_angle_cut'   : 90.0,
                   'sb_dist_opt'       : 4.0,
                   'sb_dist_cut'       : 5.0,
                   'contact_scale_opt' : 1.25,
                   'contact_scale_cut' : 1.75},
        'rd2':    {'version'           : 'rd2',
                   'level'             : 'atom',
                   'hbond_dist_opt'    : 3.5,
                   'hbond_dist_cut'    : 4.0,
                   'hbond_angle_opt'   : 60.0,
                   'hbond_angle_cut'   : 90.0,
                   'sb_dist_opt'       : 4.0,
                   'sb_dist_cut'       : 5.0,
                   'contact_scale_opt' : 1.25,
                   'contact_scale_cut' : 1.75}
      }

class Features:
    def __init__(self, root, ifp_version='rd1', shape_version='pharm_max',
                 mcss_version='mcss16', max_poses=10000):
        self.root = root
        self.ifp_version = ifp_version
        self.shape_version = shape_version
        self.mcss_version = mcss_version
        self.mcss_file = '{}/features/{}.typ'.format(os.environ['COMBINDHOME'], mcss_version)
        self.max_poses = max_poses

        self.raw = {}
        self.pair_energies = {}
        self.single_energies = {}

    def path(self, name, base=False, **kwargs):
        if name == 'root':
            return self.root

        if base:
            return '{}/{}'.format(self.root, name)
        elif name == 'rmsd':
            return '{}/docking/{}/{}_rmsd.npy'.format(self.root, basename(kwargs['pv']), basename(kwargs['pv']))
        elif name == 'gscore':
            return '{}/gscore/{}.npy'.format(self.root, basename(kwargs['pv']))
        elif name == 'ifp':
            return '{}/ifp/{}.csv'.format(self.root, basename(kwargs['pv']))
        elif name == 'ifp-pair':
            ifp1 = basename(self.path('ifp', pv=kwargs['pv1']))
            ifp2 = basename(self.path('ifp', pv=kwargs['pv2']))
            return '{}/ifp-pair/{}-{}-and-{}.npy'.format(self.root, '{}', ifp1, ifp2)
        elif name == 'shape':
            return '{}/shape/shape-{}-and-{}.npy'.format(self.root,
                basename(kwargs['pv1']), basename(kwargs['pv2']))
        elif name == 'mcss':
            return '{}/mcss/mcss-{}-and-{}.npy'.format(self.root,
                basename(kwargs['pv1']), basename(kwargs['pv2']))

    def load_features(self, features):
        self.raw = {}

        self.raw['rmsd'] = {}
        paths = self.path('rmsd', pv='*')
        paths = glob(paths)
        for path in paths:
            name = basename(path)[:-5]
            self.raw['rmsd'][name] = np.load(path)

        self.raw['gscore'] = {}
        paths = self.path('gscore', pv='*')
        paths = glob(paths)
        for path in paths:
            name = basename(path)
            self.raw['gscore'][name] = np.load(path)

        for feature in features:
            if feature in ['shape', 'mcss']:
                paths = self.path(feature, pv1='*', pv2='*')
            else:
                paths = self.path('ifp-pair', pv1='*', pv2='*').format(feature)
            paths = glob(paths)

            self.raw[feature] = {}
            for i, path in enumerate(paths):
                if not i % 100:
                    print('loading', feature, i, 'of', len(paths))
                _path = basename(path)[len(feature)+1:]
                name1, name2 = _path.split('-and-')
                self.raw[feature][(name1, name2)] = np.load(path)

    def compute_features(self, pvs, processes=1, ifp=True, shape=True, mcss=True):
        mkdir(self.path('root'))

        print('Extracting glide scores.')
        mkdir(self.path('gscore', base=True))
        for pv in pvs:
            out = self.path('gscore', pv=pv)
            if not os.path.exists(out):
                self.compute_gscore(pv, out)

        if ifp:
            print('Computing interaction fingerprints.')
            mkdir(self.path('ifp', base=True))
            unfinished = []
            for pv in pvs:
                out = self.path('ifp', pv=pv)
                if not os.path.exists(out):
                    unfinished += [(pv, out)]
            mp(self.compute_ifp, unfinished, processes)

            print('Computing interaction similarities.')
            mkdir(self.path('ifp-pair', base=True))
            unfinished = []
            for i, pv1 in enumerate(pvs):
                for pv2 in pvs[i+1:]:
                    ifp1 = self.path('ifp', pv=pv1)
                    ifp2 = self.path('ifp', pv=pv2)
                    for feature in ['hbond', 'saltbridge', 'contact']:
                        out = self.path('ifp-pair', pv1=pv1, pv2=pv2)
                        out = out.format(feature)
                        if not os.path.exists(out):
                            unfinished += [(ifp1, ifp2, feature, out)]
            mp(self.compute_ifp_pair, unfinished, processes)

        if shape:
            print('Computing shape similarities.')
            mkdir(self.path('shape', base=True))
            unfinished = []
            for i, pv1 in enumerate(pvs):
                for pv2 in pvs[i+1:]:
                    out = self.path('shape', pv1=pv1, pv2=pv2)
                    if not os.path.exists(out):
                        unfinished += [(pv1, pv2, out)]
            mp(self.compute_shape, unfinished, processes)

        if mcss:
            print('Computing mcss similarities.')
            mkdir(self.path('mcss', base=True))
            unfinished = []
            for i, pv1 in enumerate(pvs):
                for pv2 in pvs[i+1:]:
                    out = self.path('mcss', pv1=pv1, pv2=pv2)
                    if not os.path.exists(out):
                        unfinished += [(pv1, pv2, out)]
            mp(self.compute_mcss, unfinished, processes)

    def compute_gscore(self, pv, out):
        gscores = []
        with StructureReader(pv) as sts:
            next(sts)
            for st in sts:
                gscores += [st.property['r_i_docking_score']]
                if len(gscores) == self.max_poses:
                    break
        np.save(out, gscores)

    def compute_ifp(self, pv, out):
        from features.ifp import ifp
        settings = IFP[self.ifp_version]
        ifp(settings, pv, out, self.max_poses)

    def compute_ifp_pair(self, ifp1, ifp2, feature, out):
        from features.ifp_similarity import ifp_tanimoto
        tanimotos = ifp_tanimoto(ifp1, ifp2, feature)
        np.save(out, tanimotos)

    def compute_shape(self, pv1, pv2, out):
        from features.shape import shape
        sims = shape(pv2, pv1, version=self.shape_version, max_poses=self.max_poses).T
        np.save(out, sims)

    def compute_mcss(self, pv1, pv2, out):
        from features.mcss import mcss
        rmsds = mcss(pv1, pv2, self.mcss_file, self.max_poses)
        np.save(out, rmsds)

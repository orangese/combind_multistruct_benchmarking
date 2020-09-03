import os
import subprocess
import pandas as pd
import numpy as np
from glob import glob
from schrodinger.structure import StructureReader

def ifp_tanimoto(ifp1, ifp2, features):
    ifp1 = read_ifp(ifp1)
    ifp2 = read_ifp(ifp2)
    n1 = max(ifp1.index.get_level_values(level=1))+1
    n2 = max(ifp2.index.get_level_values(level=1))+1

    tanimotos = {}
    for feature in features:
        tanimotos[feature] = np.zeros((n1, n2))+0.5
        if feature not in ifp1.index.get_level_values(0):
            continue
        for i, _ifp1 in ifp1.loc[feature].groupby(['pose']):
            _ifp1 = _ifp1.droplevel(0)
            if feature not in ifp2.index.get_level_values(0):
                continue
            for j, _ifp2 in ifp2.loc[feature].groupby(['pose']):
                _ifp2 = _ifp2.droplevel(0)

                joined = _ifp1.join(_ifp2, how='outer', lsuffix='1', rsuffix='2')
                joined = joined.fillna(0)

                overlap = sum((joined['score1']*joined['score2'])**0.5)
                total = sum(joined['score1']+joined['score2'])

                tanimotos[feature][i, j] = (1 + overlap) / (2 + total - overlap)
    return tanimotos

def read_ifp(csv):
    df = pd.read_csv(csv)
    df.loc[df.label=='hbond_acceptor', 'protein_res'] = \
        [res+'acceptor' for res in df.loc[df.label=='hbond_acceptor', 'protein_res']]
    df.loc[df.label=='hbond_donor', 'protein_res'] = \
        [res+'donor' for res in df.loc[df.label=='hbond_donor', 'protein_res']]
    df.loc[df.label=='hbond_acceptor', 'label'] = 'hbond'
    df.loc[df.label=='hbond_donor', 'label'] = 'hbond'

    df = df.set_index(['label', 'pose',  'protein_res'])
    df = df.sort_index()
    return df

def basename(path):
    x = os.path.basename(path)
    x = os.path.splitext(x)[0]
    return x

def mkdir(path):
  if not os.path.exists(path):
    os.mkdir(path)

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

def get_feature(ifp, pose, feature):
    try:
        return ifp.loc[(feature, pose)]
    except:
        return ifp.iloc[:0]


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

    def filter_native(self, pv, native, thresh=2.0):
        with StructureReader(native) as sts:
            native = list(sts)
            assert len(native) == 1
            native = native[0]

        out = pv.replace('_pv.maegz', '_native_pv.maegz')
        with StructureReader(pv) as reader, StructureWriter(out) as writer:
            writer.append(next(reader))
            for st in reader:
                conf_rmsd = ConformerRmsd(native, st)
                if conf_rmsd.calculate() < thresh:
                    writer.append(st)

    def load_features(self, features):
        self.raw = {}

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
            for path in paths:
                _path = basename(path)[len(feature)+1:]
                name1, name2 = _path.split('-and-')
                self.raw[feature][(name1, name2)] = np.load(path)

    def compute_features(self, pvs, ifp=True, shape=True, mcss=True):
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
            for pv in pvs:
                out = self.path('ifp', pv=pv)
                if not os.path.exists(out):
                    self.compute_ifp(pv, out)

            print('Computing interaction similarities.')
            mkdir(self.path('ifp-pair', base=True))
            for i, pv1 in enumerate(pvs):
                for pv2 in pvs[i+1:]:
                    ifp1 = self.path('ifp', pv=pv1)
                    ifp2 = self.path('ifp', pv=pv2)
                    out = self.path('ifp-pair', pv1=pv1, pv2=pv2)
                    if not os.path.exists(out.format('saltbridge')):
                        self.compute_ifp_pair(ifp1, ifp2, out)

        if shape:
            print('Computing shape similarities.')
            mkdir(self.path('shape', base=True))
            for i, pv1 in enumerate(pvs):
                for pv2 in pvs[i+1:]:
                    out = self.path('shape', pv1=pv1, pv2=pv2)
                    if not os.path.exists(out):
                        self.compute_shape(pv1, pv2, out)

        if mcss:
            print('Computing mcss similarities.')
            mkdir(self.path('mcss', base=True))
            for i, pv1 in enumerate(pvs):
                for pv2 in pvs[i+1:]:
                    out = self.path('mcss', pv1=pv1, pv2=pv2)
                    if not os.path.exists(out):
                        self.compute_mcss(pv1, pv2, out)

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
        settings = IFP[self.ifp_version]
        settings = ' '.join('--{} {}'.format(k, v)
                            for k, v in settings.items() if k != 'version')
        cmd = '{home}/rdpython {home}/features/ifp.py {pv} {out} {max_poses} {settings}'
        cmd = cmd.format(home=os.environ['COMBINDHOME'], pv=pv, out=out,
                         max_poses=self.max_poses, settings=settings)
        subprocess.run(cmd, shell=True)

    def compute_ifp_pair(self, ifp1, ifp2, out):
        features = ['hbond', 'saltbridge', 'contact']
        tanimotos = ifp_tanimoto(ifp1, ifp2, features)
        for feature in features:
            np.save(out.format(feature), tanimotos[feature])

    def compute_shape(self, pv1, pv2, out):
        from features.shape import shape
        sims = shape(pv1, pv2, version=self.shape_version, max_poses=self.max_poses)
        np.save(out, sims)

    def compute_mcss(self, pv1, pv2, out):
        from features.mcss import mcss
        rmsds = mcss(pv1, pv2, self.mcss_file, self.max_poses)
        np.save(out, rmsds)

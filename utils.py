from multiprocessing import Pool
import os
from schrodinger.structure import StructureReader

def pv_path(root, name):
    if '_native' in name:
        name = name.replace('_native', '')
        return '{}/{}/{}_native_pv.maegz'.format(root, name, name)
    return '{}/{}/{}_pv.maegz'.format(root, name, name)

def get_pose(pv, pose):
    with StructureReader(pv) as sts:
        for _ in range(pose+1):
            next(sts)
        st = next(sts)
    return st

def basename(path):
    x = os.path.basename(path)
    x = os.path.splitext(x)[0]
    return x

def mp(function, unfinished, processes):
    if unfinished:
        with Pool(processes=processes) as pool:
            pool.starmap(function, unfinished)

def mkdir(path):
    if not os.path.exists(path):
        os.mkdir(path)

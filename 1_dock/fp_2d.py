import os
from schrodinger.structure import StructureReader
from schrodinger.application.canvas.fingerprint import CanvasFingerprintGenerator
from schrodinger.utils import log

def get_2d_fp():

    logger = log.get_output_logger('uh')

    fp_obj = CanvasFingerprintGenerator(logger, default_type='molprint2D')
    fp_obj.setPrecision(64)

    os.system('rm -rf ifp/2d_fp')
    os.system('mkdir -p ifp/2d_fp')

    for l in os.listdir('final_ligands'):
        print l
        st = StructureReader('final_ligands/{}'.format(l)).next()
        with open('ifp/2d_fp/{}.fp'.format(l.split('.')[0]), 'w') as f:
            f.write(str(fp_obj.generate(st)))

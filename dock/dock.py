from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.rmsd import ConformerRmsd
import os
import subprocess

GLIDE_ES4 = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   100
POSTDOCK_NPOSE   100
PRECISION   SP
NENHANCED_SAMPLING   4
'''

GLIDE = '''GRIDFILE  {grid}
LIGANDFILE   {ligands}
DOCKING_METHOD   confgen
POSES_PER_LIG   30
POSTDOCK_NPOSE   30
PRECISION   SP
'''

def dock(grid, ligands, root, name, enhanced, n_processes):
    infile = GLIDE_ES4 if enhanced else GLIDE
    glide_in = '{}/{}.in'.format(root, name)
    glide_pv = '{}/{}_pv.maegz'.format(root, name)
    glide_log = '{}/{}.log'.format(root, name)
    glide_cmd = 'glide -WAIT {} -HOST localhost:{}'.format(os.path.basename(glide_in), n_processes)

    if os.path.exists(glide_pv):
        return

    if os.path.exists(glide_log):
        with open(glide_log) as fp:
            logtxt = fp.read()
        phrases = ['** NO ACCEPTABLE LIGAND POSES WERE FOUND **',
                   'NO VALID POSES AFTER MINIMIZATION: SKIPPING.',
                   'No Ligand Poses were written to external file',
                   'GLIDE WARNING: Skipping refinement, etc. because rough-score step failed.']

        if any(phrase in logtxt for phrase in phrases):
            return

    if not os.path.exists(root):
        os.system('mkdir {}'.format(root))
    with open(glide_in, 'w') as fp:
        fp.write(infile.format(grid=grid, ligands=ligands))

    subprocess.run(glide_cmd, cwd=root, shell=True)

def filter_native(native, pv, out, thresh):
    with StructureReader(native) as sts:
        native = list(sts)
        assert len(native) == 1, len(native)
        native = native[0]

    
    with StructureReader(pv) as reader, StructureWriter(out) as writer:
        writer.append(next(reader))
        for st in reader:
            conf_rmsd = ConformerRmsd(native, st)
            if conf_rmsd.calculate() < thresh:
                writer.append(st)

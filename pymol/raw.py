from pymol import cmd
from glob import glob

def raw(protein):
    for prot in sorted(glob('{}/structures/aligned/*'.format(protein))):
        pdb = prot.split('/')[-1]
        load_crystal_protein(protein, pdb)
        load_crystal_pose(protein, pdb)
        
    cmd.util.cbao("prot_*")
    cmd.util.cbay("het and crystal_*")
    cmd.show('sticks', "het and crystal_*")
    cmd.hide('lines', 'element h')
    cmd.show('spheres', 'het and prot* and (crystal* expand 5)')

    cmd.show('cartoon')
    cmd.set('cartoon_oval_length', '0.5')
    cmd.set('cartoon_transparency', '0.5')
    cmd.hide('everything', 'element H and not (element N+O extend 1)')
    cmd.hide('everything', 'name H')

def load_crystal_protein(protein, ligand):
    cmd.load('{0}/structures/aligned/{1}/{1}_prot.mae'.format(protein, ligand))
    cmd.set_name('{}_prot'.format(ligand), 'prot_{}'.format(ligand))

def load_crystal_pose(protein, ligand):
    cmd.load('{0}/structures/aligned/{1}/{1}_lig.mae'.format(protein, ligand))
    cmd.set_name('{}_lig'.format(ligand), 'crystal_{}'.format(ligand))

def ambiguous(protein):
    for lig in sorted(glob('{}/structures/initial_files/*_*_lig.mae'.format(protein))):
        pdb, het, _ = lig.split('/')[-1].split('_')
        complexname = '{}-{}'.format(pdb, het)
        ligname = '{}-lig'.format(complexname)
        protname = '{}-prot'.format(complexname)
        cmd.load(lig, ligname)
        cmd.show('sticks', ligname)
        cmd.load(lig.replace('lig', 'prot'), protname)

        cmd.create(complexname, '{} or {}'.format(ligname, protname))

        cmd.delete(ligname)
        cmd.delete(protname)

def alignall(ref):
    for obj in cmd.get_names():
        cmd.align(obj, ref)

cmd.extend('raw', raw)
cmd.extend('ambiguous', ambiguous)
cmd.extend('alignall', alignall)

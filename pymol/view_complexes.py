from pymol import cmd
from glob import glob

def load_complexes(protein):
    cmd.delete('*')
    for prot in sorted(glob('{}/structures/proteins/*_prot.mae'.format(protein)))[:20]:
        protein = prot.split('/')[0]
        pdb = prot.split('/')[-1].split('_')[0]
        print pdb
        cmd.load('{}/structures/proteins/{}_prot.mae'.format(protein, pdb))
        cmd.load('{}/structures/ligands/{}_lig.mae'.format(protein, pdb))
    cmd.util.cbao("*_prot")
    cmd.util.cbay("het and *_lig")
    cmd.show('sticks', "het and *_lig")
    cmd.hide('lines', 'element h')
    cmd.show('spheres', 'het and *_prot and (*_lig expand 5)')

    cmd.show('cartoon')
    cmd.set('cartoon_oval_length', '0.5')
    cmd.set('cartoon_transparency', '0.5')
    cmd.hide('everything', 'element H and not (element N+O extend 1)')
    cmd.hide('everything', 'name H')

def show_prot():
    cmd.enable('*prot')

def show_lig():
    cmd.disable('*')
    cmd.enable('*lig')
    
cmd.extend('load_complexes', load_complexes)
cmd.extend('show_lig', show_lig)
cmd.extend('show_prot', show_prot)

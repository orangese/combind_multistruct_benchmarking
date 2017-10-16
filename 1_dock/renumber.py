import os
from schrodinger.structure import StructureReader, StructureWriter
import subprocess

script = '/share/PI/rondror/jbelk/combind/1_dock/adjust_residue_numbering.py'
schro = os.environ.get("SCHRODINGER", None)

def renumber():
    if 'renumbered_proteins' in os.listdir('.'):
        if len(os.listdir('renumbered_proteins')) == len(os.listdir('processed_proteins')):
            print 'already done'
            return
    
    os.system('rm -rf renumbered_proteins')
    os.system('mkdir renumbered_proteins')

    assert len(os.listdir('grids')) == 1
    ref_name = '{}.mae'.format(os.listdir('grids')[0])
    print ref_name, 'reference'
 
    other_names = [i for i in os.listdir('processed_proteins') if i != ref_name]
    print len(other_names), 'others'
    job_input = StructureWriter('renumbered_proteins/job_input.mae')
    ref_st = StructureReader('processed_proteins/{}'.format(ref_name)).next()

    job_input.append(ref_st)
    for n in other_names:
        other_st = StructureReader('processed_proteins/{}'.format(n)).next()
        job_input.append(other_st)
    job_input.close()

   
    os.system('{}/run {} renumbered_proteins/job_input.mae renumbered_proteins/job_output.mae -nosuper -renumber reference'.format(schro, script))

    #subprocess.call(['{}/run'.format(schro), script, 'renumbered_proteins/job_input.mae',
    #    'renumbered_proteins/job_output.mae', '-nosuper', '-renumber', 'reference'])
    
    for st in StructureReader('renumbered_proteins/job_output.mae'):
        st_out = StructureWriter('renumbered_proteins/{}.mae'.format(st._getTitle()))
        st_out.append(st)
        st_out.close()

    os.system('rm renumbered_proteins/job_input.mae renumbered_proteins/job_output.mae')
    
    print 'renumbered {} of {} structs'.format(len(os.listdir('renumbered_proteins')), len(os.listdir('processed_proteins')))
    for i in os.listdir('processed_proteins'):
        if i not in os.listdir('renumbered_proteins'): print i, 'not done'

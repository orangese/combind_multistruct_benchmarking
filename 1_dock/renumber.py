import os
import sys

from schrodinger.structure import StructureReader, StructureWriter

temp_prot = [f for f in os.listdir('.') if f.split('_')[-1] == 'template.mae' and f[:4] == 'rot-']
quer_prot = [f for f in os.listdir('.') if f.split('_')[-1] == 'query.mae' and f[:4] == 'rot-']

template = StructureReader(temp_prot[0]).next()
query = StructureReader(quer_prot[0]).next()

renumber_in = StructureWriter('renumber_in.mae')
renumber_in.append(template)
renumber_in.append(query)
renumber_in.close()

os.system('$SCHRODINGER/run adjust_residue_numbering.py -nosuper -renumber reference renumber_in.mae renumber_out.mae')

try:
    for st in StructureReader('renumber_out.mae'):
        st_wr = StructureWriter('{}_out.mae'.format(st.title))
        st_wr.append(st)
        st_wr.close()
except:
    print 'failed'

query_pdb = quer_prot[0].split('_')[0].split('-')[-1]

if not os.path.exists('{}_out.mae'.format(query_pdb)):
    print 'renumber failed', query_pdb
    os.system('cp {} {}_out.mae'.format(quer_prot[0], query_pdb))

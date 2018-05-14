import os
import sys

from schrodinger.structure import StructureReader, StructureWriter

#command1 = '$SCHRODINGER/utilities/prepwizard -WAIT -epik_pHt 2.0 -f 3 -samplewater -keepfarwat -watdist 0 {}_in.mae {}_out.mae'
#command2 = '$SCHRODINGER/utilities/prepwizard -WAIT -nopropka -watdist 0 {}_in.mae {}_out.mae'
#command3 = '$SCHRODINGER/utilities/prepwizard -WAIT {}_in.mae {}_out.mae'

command = '$SCHRODINGER/utilities/prepwizard -WAIT -watdist 0 {}_in.mae {}_out.mae'

def load_complex(pdb_id):
    prot_in = 'structures/raw_files/{}_prot.mae'.format(pdb_id)
    lig_in = 'structures/raw_files/{}_lig.mae'.format(pdb_id)
    
    prot_st = StructureReader(prot_in).next()
    
    if not os.path.exists(lig_in): 
        prot_st.title = pdb_id
        return prot_st

    lig_st = StructureReader(lig_in).next()

    assert len(lig_st.chain) == 1, pdb_id
    for c in lig_st.chain:
        c.name = 'L'

    alpha = 'ABCDEFGHIJKMNOPQRST'
    alpha_count = 0 
    for c in prot_st.chain:
        if c.name.strip() == '': continue

        c.name = alpha[alpha_count]
        alpha_count += 1
    
    merged_st = lig_st.merge(prot_st)
    merged_st.title = pdb_id
    return merged_st

def process_structs():
    os.system('mkdir -p structures/processed_files')

    all_f = sorted(os.listdir('structures/raw_files'))
    all_prot = [f.split('_')[0] for f in all_f if f.split('_')[1] == 'prot.mae']

    for pdb_id in all_prot:
        if pdb_id[0] == '.': continue

        o_dir = '{}'.format(pdb_id)

        if os.path.exists('structures/processed_files/{}/{}_out.mae'.format(o_dir, o_dir)):
            continue
      
        os.system('rm -rf structures/processed_files/{}'.format(o_dir))
        os.system('mkdir -p structures/processed_files/{}'.format(o_dir))

        merged_st = load_complex(pdb_id)
        st_wr = StructureWriter('structures/processed_files/{}/{}_in.mae'.format(o_dir, o_dir))
        st_wr.append(merged_st)
        st_wr.close()

        os.chdir('structures/processed_files/{}'.format(o_dir))
        print 'processing', o_dir
        with open('process_in.sh', 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write(command.format(o_dir, o_dir))
        os.system('sbatch -p owners -t 00:30:00 -o slurm.out process_in.sh')
        os.chdir('../../..')

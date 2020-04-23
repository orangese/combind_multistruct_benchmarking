import os
from schrodinger.structure import StructureReader

command = '$SCHRODINGER/utilities/prepwizard -WAIT -rehtreat -watdist 0 {}_in.mae {}_out.mae'
queue = 'rondror'

def load_complex(pdb):
    prot_in = 'structures/raw/{}_prot.mae'.format(pdb)
    lig_in = 'structures/raw/{}_lig.mae'.format(pdb)
    
    prot_st = next(StructureReader(prot_in))
    
    if not os.path.exists(lig_in): 
        prot_st.title = pdb
        return prot_st

    lig_st = next(StructureReader(lig_in))

    assert len(lig_st.chain) == 1, pdb
    for c in lig_st.chain:
        c.name = 'L'

    alpha = 'ABCDEFGHIJKMNOPQRST'
    alpha_count = 0
    for c in prot_st.chain:
        if c.name.strip() == '': continue

        c.name = alpha[alpha_count]
        alpha_count += 1
    
    merged_st = lig_st.merge(prot_st)
    merged_st.title = pdb
    return merged_st

def struct_process():
    if not os.path.exists('structures/raw'):
        return

    all_f = sorted(os.listdir('structures/raw'))
    all_prot = [f.split('_')[0] for f in all_f if f.split('_')[1] == 'prot.mae']

    for pdb in all_prot:
        if pdb[0] == '.': continue

        if os.path.exists('structures/processed/{}/{}_out.mae'.format(pdb, pdb)):
            continue
     
        os.system('mkdir -p structures/processed')
        os.system('rm -rf structures/processed/{}'.format(pdb))
        os.system('mkdir structures/processed/{}'.format(pdb))

        merged_st = load_complex(pdb)
        merged_st.write('structures/processed/{}/{}_in.mae'.format(pdb, pdb))

        os.chdir('structures/processed/{}'.format(pdb))
        print('processing', pdb)
        with open('process_in.sh', 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(command.format(pdb, pdb))
        os.system('sbatch -p {} -t 01:30:00 -o slurm.out process_in.sh'.format(queue))
        os.chdir('../../..')

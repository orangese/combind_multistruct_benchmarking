import os
import sys

from schrodinger.structure import StructureReader, StructureWriter

#command = '$SCHRODINGER/utilities/prepwizard -WAIT -epik_pHt 2.0 -f 3 -samplewater -keepfarwat -watdist 0 {}_in.mae {}_out.mae'
command = '$SCHRODINGER/utilities/prepwizard -WAIT -nopropka {}_in.mae {}_out.mae'

def load_complex(pdb_id):
    prot_in = 'structures/raw_files/{}_prot.mae'.format(pdb_id)
    lig_in = 'structures/raw_files/{}_lig.mae'.format(pdb_id)
    
    prot_st = StructureReader(prot_in).next()
    
    if not os.path.exists(lig_in): 
        prot_st._setTitle(pdb_id)
        return prot_st

    lig_st = StructureReader(lig_in).next()

    assert len(lig_st.chain) == 1, pdb_id
    for c in lig_st.chain:
        c._setChainName('L')

    alpha = 'ABCDEFGHIJKMNOPQRST'
    alpha_count = 0 
    for c in prot_st.chain:
        if c._getChainName().strip() == '': continue

        c._setChainName(alpha[alpha_count])
        alpha_count += 1
    
    merged_st = lig_st.merge(prot_st)
    merged_st._setTitle(pdb_id)
    return merged_st

def process_structs():
    #all_files = 'structures/processed_files'
    os.system('mkdir -p structures/processed_files')

    all_f = sorted(os.listdir('structures/raw_files'))
    all_prot = [f.split('_')[0] for f in all_f if f.split('_')[1] == 'prot.mae']

    for pdb_id in all_prot:
        if pdb_id[0] == '.': continue
        #if not os.path.exists('structures/aligned_files/{}/{}_out.mae'.format(pdb_id, pdb_id)):
        #    print 'not aligned', pdb_id
        #    continue

        if os.path.exists('structures/processed_files/{}/{}_out.mae'.format(pdb_id, pdb_id)):
            continue
        
        os.system('rm -rf structures/processed_files/{}'.format(pdb_id))
        os.system('mkdir -p structures/processed_files/{}'.format(pdb_id))
        #os.system('cp structures/aligned_files/{}/{}_out.mae {}/{}/{}_in.mae'.format(pdb_id, pdb_id, all_files, pdb_id, pdb_id))

        merged_st = load_complex(pdb_id)
        st_wr = StructureWriter('structures/processed_files/{}/{}_in.mae'.format(pdb_id, pdb_id))
        st_wr.append(merged_st)
        st_wr.close()

        os.chdir('structures/processed_files/{}'.format(pdb_id))
        print 'processing', pdb_id
        with open('process_in.sh', 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger/2017-3\n')
            f.write(command.format(pdb_id, pdb_id))
        os.system('sbatch -p rondror -t 00:30:00 -o slurm.out process_in.sh')
        os.chdir('../../..')

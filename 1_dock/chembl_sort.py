import os
import sys
from grouper import grouper

from schrodinger.structure import StructureReader, StructureWriter

from parse_chembl import load_chembl_raw, load_chembl_proc

queue = 'rondror'
group_size = 10

def get_ligands():
    os.system('mkdir -p ligands/raw_files')

    for f_name in os.listdir('structures/ligands'):
        if f_name.split('_')[1] != 'lig.mae': continue
        if os.path.exists('ligands/raw_files/{}'.format(f_name)): continue
        os.system('cp structures/ligands/{} ligands/raw_files/{}'.format(f_name, f_name))

    ligs = load_chembl_raw()
    if len(ligs) == 0: return
    written_ligs = load_chembl_proc()

    with open('chembl/chembl_info.txt', 'a') as f:
        for lig_name in sorted(ligs.keys(), key=lambda x: ligs[x].ki):
            if lig_name+'_lig' not in written_ligs:
                #if os.path.exists('ligands/raw_files/{}_lig.mae'.format(lig_name)): 
                    #os.system('rm ligands/raw_files/{}_lig.mae'.format(lig_name)) 
                    #continue
                valid_stereo = ligs[lig_name].check_stereo()
                f.write('{},{},{},{},{}\n'.format(lig_name, ligs[lig_name].target_prot, valid_stereo, 
                    ligs[lig_name].ki, ligs[lig_name].smiles))
                st_writer = StructureWriter('ligands/raw_files/{}_lig.mae'.format(lig_name))
                st_writer.append(ligs[lig_name].st)
                st_writer.close()

def proc_ligands():
    add_h = '$SCHRODINGER/utilities/prepwizard -WAIT -noepik -noprotassign -noimpref {}_in.mae {}_in_epik.mae\n' 
    epik_command = '$SCHRODINGER/epik -WAIT -ph 7.0 -pht 2.0 -imae {}_in_epik.mae -omae {}_out.mae\n'

    os.system('mkdir -p ligands/prepared_ligands')
    os.system('rm -f ligands/prepared_ligands/*.sh')
    os.system('rm -f ligands/prepared_ligands/*.out')

    all_u = [l.split('.')[0] for l in os.listdir('ligands/raw_files')]
    all_u = [l for l in all_u if not os.path.exists('ligands/prepared_ligands/{}/{}.mae'.format(l,l))]

    if len(all_u) > 0:
        print len(all_u), 'unfinished ligands'

    unfinished = []
    for l in all_u:
        prepped = 'ligands/prepared_ligands/{}/{}_out.mae'.format(l,l)
        final = 'ligands/prepared_ligands/{}/{}.mae'.format(l,l)
        if os.path.exists(prepped): 
            if not os.path.exists(scaff):
                try:
                    st = StructureReader(prepped).next() # first protonation state
                    st.title = l
                    stwr = StructureWriter(final)
                    stwr.append(st)
                    stwr.close()
                except Exception as e:
                    print l
                    print e
                    break
        else:
            unfinished.append(l)

    if len(unfinished) > 0:
        print len(unfinished), 'unprocessed ligands'
    
    for i, ligs in enumerate(grouper(group_size, unfinished)):
        with open('ligands/prepared_ligands/batch-{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            for name in ligs:
                if name is None: continue

                os.system('rm -rf ligands/prepared_ligands/{}'.format(name))
                os.system('mkdir ligands/prepared_ligands/{}'.format(name))
                os.system('cp ligands/raw_files/{}.mae ligands/prepared_ligands/{}/{}_in.mae'.format(name, name, name))

                f.write('cd {}\n'.format(name))
                f.write('sh process_in.sh > process.out\n')
                f.write('cd ..\n')

                with open('ligands/prepared_ligands/{}/process_in.sh'.format(name), 'w') as f2:
                    f2.write('#!/bin/bash\nmodule load schrodinger\n')
                    f2.write(add_h.format(name, name))
                    f2.write(epik_command.format(name, name))

        os.chdir('ligands/prepared_ligands')
        os.system('sbatch -p {} -t 1:00:00 batch-{}.sh'.format(queue,i))
        os.chdir('../..')




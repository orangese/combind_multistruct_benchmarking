import os
import sys

from schrodinger.structure import SmilesStructure, StructureReader, StructureWriter

from parse_chembl import load_chembl_raw, load_chembl_proc, desalt, load_drugs

def get_drugs():
    os.system('mkdir -p ligands/drugs')

    drugs = load_drugs() 
    for c_id, smi in drugs.items():
        if not os.path.exists('ligands/raw_files/{}_lig.mae'.format(c_id)):
            smi_st = desalt(SmilesStructure(smi).get3dStructure(False))
            smi_st._setTitle('{}_lig'.format(c_id))
            st_wr = StructureWriter('ligands/raw_files/{}_lig.mae'.format(c_id))
            st_wr.append(smi_st)
            st_wr.close()

def get_ligands():
    os.system('mkdir -p ligands/raw_files')

    ligs = load_chembl_raw()
    written_ligs = load_chembl_proc()

    for f_name in os.listdir('structures/raw_files'):
        if f_name.split('_')[1] != 'lig.mae': continue
        if os.path.exists('ligands/raw_files/{}'.format(f_name)): continue
        os.system('cp structures/raw_files/{} ligands/raw_files/{}'.format(f_name, f_name))

    with open('chembl/chembl_info.txt', 'a') as f:
        for lig_name in sorted(ligs.keys(), key=lambda x: ligs[x].ki):
            if lig_name not in written_ligs:
                if os.path.exists('ligands/raw_files/{}_lig.mae'.format(lig_name)): 
                    os.system('rm ligands/raw_files/{}_lig.mae'.format(lig_name)) 
                    #continue
                valid_stereo = ligs[lig_name].check_stereo()
                f.write('{},{},{},{},{}\n'.format(lig_name, ligs[lig_name].target_prot, valid_stereo, ligs[lig_name].ki, ligs[lig_name].smiles))
                st_writer = StructureWriter('ligands/raw_files/{}_lig.mae'.format(lig_name))
                st_writer.append(ligs[lig_name].st)
                st_writer.close()

def proc_ligands():
    add_h = '$SCHRODINGER/utilities/prepwizard -WAIT -noepik -noprotassign -noimpref {}_in.mae {}_in_epik.mae\n' 
    epik_command = '$SCHRODINGER/epik -WAIT -ph 7.0 -pht 2.0 -imae {}_in_epik.mae -omae {}_out.mae\n'

    group_size = 15
    group_index = 0

    os.system('mkdir -p ligands/prepared_ligands')

    os.system('rm -f ligands/prepared_ligands/*.sh')
    os.system('rm -f ligands/prepared_ligands/*.out')

    for l in os.listdir('ligands/raw_files'):
        name = l.split('.')[0]
        if os.path.exists('ligands/prepared_ligands/{}/{}_out.mae'.format(name, name)):
            continue

        os.system('rm -rf ligands/prepared_ligands/{}'.format(name))
        os.system('mkdir ligands/prepared_ligands/{}'.format(name))

        os.system('cp ligands/raw_files/{} ligands/prepared_ligands/{}/{}_in.mae'.format(l, name, name))

        with open('ligands/prepared_ligands/batch-{}.sh'.format(group_index/group_size),'a') as f:
            if group_index % group_size == 0:
                f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write('cd {}\n'.format(name))
            f.write('sh process_in.sh > slurm.out\n')
            f.write('cd ..\n')

        with open('ligands/prepared_ligands/{}/process_in.sh'.format(name), 'w') as f:
            f.write('#!/bin/bash\nmodule load schrodinger\n')
            f.write(add_h.format(name, name))
            f.write(epik_command.format(name, name))

        if group_index % group_size == group_size - 1:# == group_size:

            os.chdir('ligands/prepared_ligands')
            os.system('sbatch -p owners -t 1:00:00 batch-{}.sh'.format(group_index/group_size))
            os.chdir('../..')
        group_index += 1

    if group_index % group_size != 0 and os.path.exists('ligands/prepared_ligands/batch-{}.sh'.format(group_index/group_size)):
        os.chdir('ligands/prepared_ligands')
        os.system('sbatch -p owners -t 1:00:00 batch-{}.sh'.format(group_index/group_size))
        os.chdir('../..')



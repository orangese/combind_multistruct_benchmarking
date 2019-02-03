import os
import sys
from grouper import grouper

from schrodinger.structure import StructureReader, StructureWriter

from dock.parse_chembl import load_chembl_raw, load_chembl_proc, load_dude_raw

queue = 'owners'
group_size = 10

#######################################################################################
# Generate unprocessed ligand MAE files

def copy_pdb_ligs():
    """
    Copies pdb ligands into ligands/raw_files.
    If ligand is already copied, does nothing.
    Called from datadir root.
    """
    for f_name in os.listdir('structures/ligands'):
        if f_name.split('_')[1] != 'lig.mae': continue
        if os.path.exists('ligands/raw_files/{}'.format(f_name)): continue
        os.system('cp structures/ligands/{} ligands/raw_files/{}'.format(f_name, f_name))

def copy_pdb_ligs_from_read_root(read_root, write_root):
    """
    Copies pdb ligands into ligands/raw_files.
    If ligand is already copied, does nothing.
    Called from datadir root.
    """
    print("Copying the following PDB ligand files:")
    for f_name in os.listdir('{}/structures/ligands'.format(read_root)):
        if f_name.split('_')[1] != 'lig.mae': continue
        if os.path.exists('{}/ligands/raw_files/{}'.format(write_root,f_name)): continue
        print("* {}".format(f_name))
        os.system('cp {}/structures/ligands/{} {}/ligands/raw_files/{}'.format(read_root, f_name, write_root, f_name))

def write_unprocessed_chembl_ligs(ligs, written_ligs):
    """
    Writes unprepped MAE files for chembl ligands and
    writes important properties to chembl_info.txt.
    """
    print("Writing the following ligands to chembl/chembl_info.txt and adding structure files to ligands/raw_files/:")
    with open('chembl/chembl_info.txt', 'a+') as f:
        for i, lig_name in enumerate(sorted(ligs.keys(), key=lambda x: ligs[x].ki)):
            if lig_name+'_lig' not in written_ligs:
                print("* decoy {}: {}".format(i, lig_name))
                valid_stereo = ligs[lig_name].check_stereo()
                f.write('{},{},{},{},{}\n'.format(lig_name, ligs[lig_name].target_prot,
                                                  valid_stereo, ligs[lig_name].ki,
                                                  ligs[lig_name].smiles))
                st_writer = StructureWriter('ligands/raw_files/{}_lig.mae'.format(lig_name))
                st_writer.append(ligs[lig_name].st)
                st_writer.close()

def get_ligands(read_root, write_root):
    """
    Generates unprocessed ligand mae files for pdb ligands(by copying from
    structures/ligands) and chembl ligands (by reading from the chembl/*.xls files).
    Also adds entries in chembl_info.txt for chembl ligands that are added.
    
    * Terminating execution will corrupt chembl_info.txt *
    
    Called from datadir root.
    """
    os.system('mkdir -p ligands/raw_files')
    os.system('mkdir -p dude/')
    os.system('mkdir -p chembl/')
    # copy_pdb_ligs()
    copy_pdb_ligs_from_read_root(read_root, write_root)
    # ligs = load_chembl_raw()
    ligs = load_dude_raw()
    if len(ligs) == 0: return
    written_ligs = load_chembl_proc()
    write_unprocessed_chembl_ligs(ligs, written_ligs)

#######################################################################################
# Process ligand files

def resolve_protonation_state(all_u):
    """
    Resolve multiple protonation states for each ligand by taking the
    first entry in the file (corresponding to the best scoring species).
    """
    count = 0
    print("Resolving protonation states for the following ligands:")
    for l in all_u:
        prepped = 'ligands/prepared_ligands/{}/{}_out.mae'.format(l,l)
        final = 'ligands/prepared_ligands/{}/{}.mae'.format(l,l)
        if os.path.exists(prepped):
            if not os.path.exists(final):
                try:
                    st = next(StructureReader(prepped)) # first protonation state
                    st.title = l
                    stwr = StructureWriter(final)
                    stwr.append(st)
                    stwr.close()
                    print("* {}".format(l))
                    count += 1
                except Exception as e:
                    print(l)
                    print(e)
                    break
    if count: print("Resolved protonation state for {} ligands.".format(count))


def run_ligand_processing(unfinished):
    """
    Run prepwizard and epik on ligands in list unfinished.
    """
    add_h = '$SCHRODINGER/utilities/prepwizard -WAIT -rehtreat -noepik -noprotassign -noimpref {}_in.mae {}_in_epik.mae\n' 
    epik_command = '$SCHRODINGER/epik -WAIT -ph 7.0 -pht 2.0 -imae {}_in_epik.mae -omae {}_out.mae\n'

    print("Running ligand processing for the following ligands:")
    for i, ligs in enumerate(grouper(group_size, unfinished)):
        for ligand in ligs:
            print("* {}".format(ligand))
        with open('ligands/prepared_ligands/batch-{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            for name in ligs:
                if name is None: continue

                os.system('rm -rf ligands/prepared_ligands/{}'.format(name))
                os.system('mkdir ligands/prepared_ligands/{}'.format(name))
                os.system('cp ligands/raw_files/{}.mae ligands/prepared_ligands/{}/{}_in.mae'.format(name, name, name))

                f.write('cd {}\n'.format(name))
                f.write('sh process_in.sh > process.out\n')
                f.write('cd ..\n')

                with open('ligands/prepared_ligands/{}/process_in.sh'.format(name), 'w') as f2:
                    f2.write('#!/bin/bash\n')
                    f2.write(add_h.format(name, name))
                    f2.write(epik_command.format(name, name))

        os.chdir('ligands/prepared_ligands')
        os.system('sbatch -p {} -t 1:00:00 batch-{}.sh'.format(queue,i))
        os.chdir('../..')

def proc_ligands():
    """
    Runs prepwizard and epik on ligands in raw_files.
    Writes output to prepared_ligands.
    """
    os.system('mkdir -p ligands/prepared_ligands')
    os.system('rm -f ligands/prepared_ligands/*.sh')
    os.system('rm -f ligands/prepared_ligands/*.out')

    # All ligands that have not been fully prepared.
    all_u = [l.split('.')[0] for l in os.listdir('ligands/raw_files')]
    all_u = [l for l in all_u
             if not os.path.exists('ligands/prepared_ligands/{}/{}.mae'.format(l,l))]
    print(len(all_u), 'unfinished ligands')

    if len(all_u) > 0:
        resolve_protonation_state(all_u)

    unfinished = []
    for l in all_u:
        prepped = 'ligands/prepared_ligands/{}/{}_out.mae'.format(l,l)
        if not os.path.exists(prepped):
            unfinished.append(l)
    print(len(unfinished), 'unprocessed ligands')
    
    if len(unfinished) > 0:
        run_ligand_processing(unfinished)

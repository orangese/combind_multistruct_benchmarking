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
    """ Copies pre-existing pdb ligand .mae files from structures/ligands/
    into ligands/raw_files
    """
    for f_name in os.listdir('structures/ligands'):
        if f_name.split('_')[1] != 'lig.mae': continue
        if os.path.exists('ligands/raw_files/{}'.format(f_name)): continue
        os.system('cp structures/ligands/{} ligands/raw_files/{}'.format(f_name, f_name))

def write_unprocessed_chembl_ligs(ligs, written_ligs):
    """ For each ligand not in chemb_info.txt (as determined by load_chembl_proc()),
    appends a line describing the ligand to chembl_info.txt. In addition, adds a raw .mae
    file for the ligand to the ligands/raw_files/ directory. The latter is done using
    Schrodinger's StructureWriter class.
    """
    print("Writing the following unprocessed ligands to chembl/chembl_info.txt:")
    with open('chembl/chembl_info.txt', 'a') as f:
        for lig_name in sorted(ligs.keys(), key=lambda x: ligs[x].ki):
            if lig_name+'_lig' not in written_ligs:
                print("* {}".format(lig_name))
                valid_stereo = ligs[lig_name].check_stereo()
                # Write a line of properties to chembl_info.txt
                f.write('{},{},{},{},{}\n'.format(lig_name, ligs[lig_name].target_prot,
                                                  valid_stereo, ligs[lig_name].ki,
                                                  ligs[lig_name].smiles))
                # Add a .mae file to ligands/raw_files/
                st_writer = StructureWriter('ligands/raw_files/{}_lig.mae'.format(lig_name))
                st_writer.append(ligs[lig_name].st)
                st_writer.close()

def get_ligands():
    """ Generates unprocessed ligand mae files for pdb ligands (by copying from
    structures/ligands) and chembl ligands (by reading from the chembl/chembl_info.txt file).
    Also adds entries in chembl_info.txt for any new chembl ligands that are loaded.
    
    * Terminating execution will corrupt chembl_info.txt *
    """
    print("Getting unprocessed ligands")
    print("---------------------------")

    os.system('mkdir -p ligands/raw_files')

    print("Copying PDB ligands from structures/ligands/ to ligands/raw_files/")

    # Copy PDB ligands from structures/ligands/ to ligands/raw_files/
    copy_pdb_ligs()

    print("Loading chembl ligands from chembl/*xls files")

    # Creates CHEMBL objects (see parse_chembl.py) for ligands in chembl/*.xls files
    ligs = load_chembl_raw() 

    print("Loading DUD-E ligands from DUD-E/*.ism files (all actives and up to 100 decoys are loaded")

    # Creates CHEMBL objects (see parse_chembl.py) for dude ligands
    ligs = load_dude_raw(ligs)
    if len(ligs) == 0: return

    # If a CHEMBL/DUDE ligand doesn't already appear in chembl_info.txt, add it, and write its
    # .mae structure file to ligands/raw_files/
    written_ligs = load_chembl_proc()
    write_unprocessed_chembl_ligs(ligs, written_ligs)

#######################################################################################
# Process ligand files

def resolve_protonation_state(all_u):
    """ If there exists an prepped, non-final .mae file for the ligand, then resolve
    the protonation process by taking the first entry in this file (corresponding to the
    best scoring species)
    """
    print("Following ligands have been finished:")
    count = 0
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
                    count += 1
                    print("* {}".format(l))
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

    print("Following ligands have been submitted for prepwizard+epik processing:")

    for i, ligs in enumerate(grouper(group_size, unfinished)):
        with open('ligands/prepared_ligands/batch-{}.sh'.format(i),'w') as f:
            f.write('#!/bin/bash\n')
            for name in ligs:
                if name is None: continue
                print("* {}".format(name))

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
    print("Processing ligands")
    print("------------------")

    os.system('mkdir -p ligands/prepared_ligands')
    os.system('rm -f ligands/prepared_ligands/*.sh')
    os.system('rm -f ligands/prepared_ligands/*.out')

    # all_u is a list of ligands appearing in ligands/raw_files/ without a corresponding
    # entry in ligands/prepared_ligands/
    all_u = [l.split('.')[0] for l in os.listdir('ligands/raw_files')]
    all_u = [l for l in all_u
             if not os.path.exists('ligands/prepared_ligands/{}/{}.mae'.format(l,l))]
    if len(all_u) > 0:
        print('{} ligands in ligands/raw_files/ are unfinished'.format(len(all_u)))

    # For each ligand, if the ligand is prepped but not finalized, finalize it now by grabbing
    # the highest scoring protonation state
    resolve_protonation_state(all_u)

    # Get all ligands that have not been prepped/finalized
    unfinished = []
    for l in all_u:
        prepped = 'ligands/prepared_ligands/{}/{}_out.mae'.format(l,l)
        if not os.path.exists(prepped):
            unfinished.append(l)
    if len(unfinished) > 0:
        print('{} ligands in ligands/raw_files/ are unprocessed'.format(len(unfinished)))
    
    run_ligand_processing(unfinished)

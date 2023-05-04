import os
from schrodinger.structure import StructureReader
import shutil

#hack to not include ligands here!
def should_include_lig(struct):
    return True #False
    # if "_" in struct:
    #     print('modeled_struct!')
    #     return True
    # return True

def split_complex(st, pdb_id):
    os.system('mkdir -p structures/proteins structures/ligands')
    lig_path = 'structures/ligands/{}_lig.mae'.format(pdb_id)
    prot_path = 'structures/proteins/{}_prot.mae'.format(pdb_id)
    print('split?', prot_path, lig_path)
    #need to copy from raw, because wasn't minimized! I'm also copyin unrotated structure here!

    # if not os.path.exists(lig_path) and len([a.index for a in st.atom if a.chain == 'L']) > 0:
    print("split_complex in struct_sort", pdb_id, should_include_lig(pdb_id))
    if not should_include_lig(pdb_id):
        # places ligand from raw/[]_lig.mae into structures/ligands/[]_lig.mae
        # this should happen if you DIDN'T have a ligand in the original raw/[]_prot.mae file,
        # and you didn't minimize around the nonexistent ligand in that file
        # this should also happen if you minimized around a DIFFERENT LIGAND during
        # the struct_process schrodinger step, and want to dock to raw/{}_lig.mae instead
        # BOTTOM LINE: USE THIS BRANCH IF YOU HAVE A LIGAND IN raw/[]_lig.mae 
        # THAT YOU WANT TO DOCK TO, AND IT'S DISTINCT FROM THE LIGAND ATTACHED
        # TO processed/[]_out.mae
        if os.path.exists('structures/raw/{}_lig.mae'.format(pdb_id)):
            shutil.copy( 'structures/raw/{}_lig.mae'.format(pdb_id), 'structures/ligands/{}_lig.mae'.format(pdb_id))
    else:
        # this should happen if the structure already has a ligand THAT YOU WANT
        # TO DOCK TO (i.e., should NOT happen for glosa and other lig exps)
        # for those other cases, make sure that structures/raw/[]_lig has the ligand you want
        # to dock to!
        # BOTTOM LINE: use if you didn't minimize during structprep, and the ligand
        # currently in processed/[]_out.mae is the one you want to dock
        lig_st = st.extract([a.index for a in st.atom if a.chain == 'L'])
        lig_st.title = '{}_lig'.format(pdb_id)
        lig_st.write(lig_path)
    
    # if not os.path.exists(prot_path):
    prot_st = st.extract([a.index for a in st.atom if a.chain != 'L'])
    prot_st.title = '{}_prot'.format(pdb_id)
    prot_st.write(prot_path)

#version used rotation, but I don't want that!
# def struct_sort(structs):
#     for struct in structs:
#         opt_complex = 'structures/aligned/{}/rot-{}_query.mae'.format(struct, struct)

#         if os.path.exists(opt_complex):
#             comp_st = next(StructureReader(opt_complex))
#             split_complex(comp_st, struct)

def struct_sort(structs):
    for struct in structs:
        #using unrotated query here. makes it easier to align everythin later. 
        # opt_complex = 'structures/aligned/{}/{}_query.mae'.format(struct, struct)
        opt_complex = 'structures/processed/{}/{}_out.mae'.format(struct, struct)
        if os.path.exists(opt_complex):
            print('help')
            comp_st = next(StructureReader(opt_complex))
            split_complex(comp_st, struct)

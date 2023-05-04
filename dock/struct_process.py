import shutil
import os
from schrodinger.structure import StructureReader
from subprocess import run

#hacky switch for min and no min - sorry joe
command_no_min = '$SCHRODINGER/utilities/prepwizard -WAIT -noimpref -rehtreat -watdist 0 {}_in.mae {}_out.mae'
command_min = '$SCHRODINGER/utilities/prepwizard -WAIT -rehtreat -watdist 0 {}_in.mae {}_out.mae'
command = command_min

#hack to not include ligands here!
def should_include_lig(struct):
    return True #False
    # if "_" in struct:
    #     print('modeled_struct!')
    #     return False
    # return True

def load_complex(prot_in, lig_in, struct):

    prot_st = next(StructureReader(prot_in))
    prot_st.title = struct

    # if raw/[]_lig.mae doesn't exist, then don't include it!
    if not os.path.exists(lig_in): 
        print(f"{lig_in=} doesn't exist, exiting with struct_to_use=None")
        prot_st.title = struct
        return prot_st, None

    lig_st = next(StructureReader(lig_in))

    assert len(lig_st.chain) == 1, struct
    for c in lig_st.chain:
        c.name = 'L'

    alpha = 'ABCDEFGHIJKMNOPQRST'
    alpha_count = 0
    for c in prot_st.chain:
        if c.name.strip() == '': continue

        c.name = alpha[alpha_count]
        alpha_count += 1
    struct_to_use = None

    # if not should_include_lig, then do a dummy merge (assume
    # the correct ligand is already in the protein, or you don't want to minimize)
    # note: this is the same behavior as if the raw/[]_lig.mae
    # file just never existed in the first place
    if not should_include_lig(struct):
        # this should happen iff
        # 1) there is NO LIGAND in the prot_in file and you want to keep it that way,
        # i.e., you don't want to run minimization around the ligand
        # 2) there is LIGAND in prot_in file, and you want to minimize around it
        print(f"using struct_to_use as struct from {prot_in}, title={struct}")
        prot_st.title = struct
        struct_to_use = prot_st.copy()

    merged_st = lig_st.merge(prot_st.copy())
    merged_st.title = struct

    return (merged_st, struct_to_use)

def struct_process(structs,
                   protein_in='structures/raw/{pdb}_prot.mae',
                   ligand_in='structures/raw/{pdb}_lig.mae',
                   processed_in_to_use = 'structures/processed/{pdb}/{pdb}_in.mae',
                   processed_out='structures/processed/{pdb}/{pdb}_out.mae',
                   processed_sh='structures/processed/{pdb}/process.sh'):
    processed_in_merged='structures/processed/{pdb}/{pdb}_merged.mae'
    for struct in structs:
        _protein_in = protein_in.format(pdb=struct)
        _ligand_in = ligand_in.format(pdb=struct)
        _processed_in_to_use = processed_in_to_use.format(pdb=struct)
        _processed_in_merged = processed_in_merged.format(pdb=struct)
        _processed_out = processed_out.format(pdb=struct)
        _processed_sh = processed_sh.format(pdb=struct)
        _workdir = os.path.dirname(_processed_sh)

        # if os.path.exists(_processed_out):
            # continue

        print('processing', struct)

        os.system('mkdir -p {}'.format(os.path.dirname(_workdir)))
        os.system('rm -rf {}'.format(_workdir))
        os.system('mkdir {}'.format(_workdir))

        (merged_st, struct_to_use)= load_complex(_protein_in, _ligand_in, struct)
        # merged_st has both the protein and ligan, saved in processed/{}_in.mae
        
        merged_st.write(_processed_in_merged)
        
        if  struct_to_use is None:
            print(f"writing merged_st to {_processed_in_to_use}")
            merged_st.write(_processed_in_to_use)
        else:
            print(f"writing struct_to_use to {_processed_in_to_use}")
            struct_to_use.write(_processed_in_to_use)

        print(f"struct_process, {_protein_in=}, {_ligand_in=}, {command.format(struct, struct)=}")
        with open('{}/process_in.sh'.format(_workdir), 'w') as f:
            f.write('#!/bin/bash\n')
            f.write(command.format(struct, struct))
        run('sh process_in.sh', shell=True, cwd=_workdir)

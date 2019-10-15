import os
import sys
from schrodinger.structure import StructureReader, StructureWriter
from dock.renumber import renumber

out_dir = 'structures/aligned_files'
queue = 'owners'

def align_successful(out_dir, struct, verbose=False):
    if struct in ['5IRX','5IS0','3J5Q']: 
        in_f = 'structures/processed_files/{}/{}_out.mae'.format(struct, struct)
        out_f = '{}/{}/{}_out.mae'.format(out_dir, struct, struct)
        if not os.path.exists(out_f):
            os.system('mkdir -p {}/{}'.format(out_dir, struct))
            os.system('cp {} {}'.format(in_f, out_f))
        return True # TRPV1 was manually aligned because 3J5Q has no ligand

    if not os.path.exists('{}/{}/rot-{}_query.mae'.format(out_dir, struct, struct)):
        return False
    
    if os.path.exists('{}/{}/{}_template.mae'.format(out_dir, struct, struct)):
        return True # query = template so we don't need to check alignment

    with open('{}/{}/align.out'.format(out_dir, struct), 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if len(tmp) > 0 and tmp[0] == 'Alignment':
                if float(tmp[2]) > 0.15:
                    print('-- Alignment warning!', struct, float(tmp[2]))
                    return False
                return True
            if verbose and len(tmp) > 0 and tmp[0] == 'RMSD:':
                print('{} to {} Alignment RMSD: {}'.format(struct, template, tmp[1]))
        else:
            print('alignment failure', struct)
            return False

def align_structs(verbose=False):

    os.system('mkdir -p {}'.format(out_dir))

    all_prot = sorted([p for p in os.listdir('structures/processed_files') if p[0] != '.'])
    template = all_prot[0]
   
    if template == '2R4R': template = '2RH1'
 
    if os.path.exists(out_dir) and len(os.listdir(out_dir)) > 0:
        for f in os.listdir('{}/{}'.format(out_dir, os.listdir(out_dir)[0])):
            if f.split('_')[-1] == 'template.mae':
                prot = f.split('_')[0]
                template = prot.split('-')[-1]
                break

    template_path = 'structures/processed_files/{}/{}_out.mae'.format(template, template)
    if not os.path.exists(template_path):
        print('template not processed', template_path)
        return

    for struct in all_prot:
        query_path = 'structures/processed_files/{}/{}_out.mae'.format(struct, struct)
        if not os.path.exists(query_path):
            continue
        if align_successful(out_dir, struct, verbose):

            if not os.path.exists('{}/{}/{}_out.mae'.format(out_dir, struct, struct)):
                print('renumber', struct)
                os.chdir('{}/{}'.format(out_dir, struct))
                renumber()
                os.chdir('../../..')
            
            continue

        os.system('rm -rf {}/{}'.format(out_dir, struct))
        os.system('mkdir -p {}/{}'.format(out_dir, struct))

        os.system('cp {} {}/{}/{}_template.mae'.format(template_path, out_dir, struct, template))
        os.system('cp {} {}/{}/{}_query.mae'.format(query_path, out_dir, struct, struct))
        
        os.chdir('{}/{}'.format(out_dir, struct))
        print('align', struct)
                
        with open('align_in.sh', 'w') as f:
            f.write('#!/bin/bash\n'
                    '$SCHRODINGER/utilities/structalign \\\n'
                    '  -asl        "(not chain. L and not atom.element H) '
                                'and (fillres within 15.0 chain. L)" \\\n'
                    '  -asl_mobile "(not chain. L and not atom.element H) '
                                'and (fillres within 15.0 chain. L)" \\\n'
                    '  {}_template.mae {}_query.mae\n'.format(template, struct))
        os.system('sbatch -p {} -t 00:10:00 -o align.out align_in.sh'.format(queue))
        os.chdir('../../..')

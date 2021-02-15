import os

out_dir = 'structures/aligned'

def align_successful(out_dir, struct):

    if not os.path.exists('{}/{}/rot-{}_query.mae'.format(out_dir, struct, struct)):
        return False
    
    if os.path.exists('{}/{}/{}_template.mae'.format(out_dir, struct, struct)):
        return True # query = template so we don't need to check alignment

    with open('{}/{}/align.out'.format(out_dir, struct), 'r') as f:
        for line in f:
            tmp = line.strip().split()
            if len(tmp) > 0 and tmp[0] == 'Alignment':
                if float(tmp[2]) > 0.4:
                    print('-- Alignment warning!', struct, float(tmp[2]))
                    return False
                return True
        else:
            print('alignment failure', struct)
            return False

def struct_align(template=None):
    if not os.path.exists('structures/processed'):
        return

    all_prot = sorted([p for p in os.listdir('structures/processed') if p[0] != '.'])
    if not len(all_prot):
        return
    if template is None:
        template = all_prot[0]
 
    if os.path.exists(out_dir) and len(os.listdir(out_dir)) > 0:
        for f in os.listdir('{}/{}'.format(out_dir, os.listdir(out_dir)[0])):
            if f.split('_')[-1] == 'template.mae':
                prot = f.split('_')[0]
                template = prot.split('-')[-1]
                break

    template_path = 'structures/processed/{}/{}_out.mae'.format(template, template)
    if not os.path.exists(template_path):
        print('template not processed', template_path)
        return

    for struct in all_prot:
        query_path = 'structures/processed/{}/{}_out.mae'.format(struct, struct)
        if not os.path.exists(query_path) or align_successful(out_dir, struct):
            continue

        os.system('mkdir -p {}'.format(out_dir))
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
        os.system('sh align_in.sh > align.out')
        os.chdir('../../..')

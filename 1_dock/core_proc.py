import os
import sys

def load_matches(data_dir=None):
    if data_dir is not None:
        data_dir = '{}/docking/core'.format(data_dir)
    else:
        data_dir = 'docking/core'
    matches = {}
    for ss in [f for f in os.listdir(data_dir) if f[-3:] == 'mae']:
        matches[ss] = {}
        cache = '{}/{}_matches.txt'.format(data_dir, ss.split('.')[0])
        if os.path.exists(cache):
            with open(cache) as f:
                for line in f:
                    lig, match = line.strip().split(',')
                    if match == 'True':
                        matches[ss][lig] = True
                    else:
                        matches[ss][lig] = False
    return matches

def core_proc():

    from schrodinger.structure import StructureReader, StructureWriter
    from schrodinger.structutils.analyze import evaluate_smarts_canvas, generate_smarts_canvas

    matches = load_matches()
    for ss in matches:
        cache = 'docking/core/{}_matches.txt'.format(ss.split('.')[0])

        smarts = generate_smarts_canvas(StructureReader('docking/core/{}'.format(ss)).next())    

        for i,lig in enumerate(os.listdir('ligands/unique')):
            lig = lig.split('.')[0]

            if lig not in matches[ss]:
                st = StructureReader('ligands/unique/{}.mae'.format(lig)).next()
                match = len(evaluate_smarts_canvas(st, smarts)) > 0
                with open(cache,'a') as f:
                    f.write('{},{}\n'.format(lig, match))

            #elif matches[ss][lig] and ss == 'HALPss2.mae':
            #    st_list = [st for st in StructureReader('ligands/prepared_ligands/{}/{}_out.mae'.format(lig, lig))]
            #    if len(st_list) > 1:
            #        print lig, len(st_list), [s.formal_charge for s in st_list]

        print ss, smarts, len([i[0] for i in matches[ss].items() if i[1]])
    

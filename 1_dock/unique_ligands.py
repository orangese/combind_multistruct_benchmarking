import os

from parse_chembl import load_chembl_proc

from schrodinger.structure import StructureReader, StructureWriter

def filter_duplicates():
    os.system('mkdir -p ligands/unique ligands/duplicate')
    
    ligs = load_chembl_proc()

    duplicates = {} 
    lig_st = {}
    new_ligs_seen = 0
    for lig_id in sorted([l.split('.')[0] for l in os.listdir('ligands/prepared_ligands')]):
        lig_path = 'ligands/prepared_ligands/{}/{}_out.mae'.format(lig_id, lig_id)
        if os.path.exists('ligands/unique/{}.mae'.format(lig_id)): continue
        if os.path.exists('ligands/duplicate/{}.mae'.format(lig_id)): continue
        if not os.path.exists(lig_path): continue

        if len(duplicates.keys()) == 0:
            for u_lig in [l.split('_')[0] for l in os.listdir('ligands/unique')]:
                ligpath = 'ligands/unique/{}.mae'.format(u_lig)
                try:
                    lig_st[u_lig] = StructureReader(ligpath).next()
                    duplicates[u_lig] = [u_lig]
                except:
                    os.system('rm -f {}'.format(ligpath))

        new_ligs_seen += 1
       
        try: 
            lig_st[lig_id] = [st for st in StructureReader(lig_path)][0] # take first epik state
        except:
            os.system('rm -rf ligands/prepared_ligands/{}'.format(lig_id))
            print 'check lig processing', lig_path
            continue

        for lig_id2 in duplicates:
            if lig_st[lig_id].isEquivalent(lig_st[lig_id2], True):
                duplicates[lig_id2].append(lig_id)
                break
        else:
            duplicates[lig_id] = [lig_id]

    if new_ligs_seen > 0: 
        print 'found {} new ligands'.format(new_ligs_seen)
        os.system('rm -f chembl/duplicates.txt')    

    for lig_id, dups in duplicates.items():
        pdb_ligs = sorted([l for l in dups if l[:6] != 'CHEMBL'])
        chembl_ligs = sorted([l for l in dups if l[:6] == 'CHEMBL'],key=lambda x: ligs[x].ki)

        keep = []
        if len(pdb_ligs) > 0:
            keep.append(pdb_ligs[0]) # alphabetically first non-chembl ligand

        if len(chembl_ligs) > 0:
            stereo = [l for l in chembl_ligs if ligs[l].valid_stereo]
            if len(stereo) > 0:
                keep.append(stereo[0]) # best ki ligand with specified stereochemistry
            else:
                keep.append(chembl_ligs[0]) # best ki ligand

        if len(keep) > 1:
            with open('chembl/duplicates.txt','a') as f:
                f.write('{},{}\n'.format(keep[0],keep[1]))
        
        for l in dups:
            path_write = 'ligands/duplicate/{}.mae'.format(l)
            path_delete = 'ligands/unique/{}.mae'.format(l)
            if l in keep:
                path_write, path_delete = path_delete, path_write
            os.system('rm -f {}'.format(path_delete))
            if os.path.exists(path_write): continue
            st_wr = StructureWriter(path_write)
            st_wr.append(lig_st[l])
            st_wr.close()




import os

from parse_chembl import load_chembl_proc, load_drugs

from schrodinger.structure import StructureReader, StructureWriter

def filter_duplicates():
    os.system('mkdir -p ligands/unique ligands/duplicate')
    
    ligs = load_chembl_proc()
    drugs = load_drugs()

    duplicates = {} 
    lig_st = {}
    new_ligs_seen = 0
    for lig_id in sorted([l.split('_')[0] for l in os.listdir('ligands/prepared_ligands')]):
        lig_path = 'ligands/prepared_ligands/{}_lig/{}_lig_out.mae'.format(lig_id, lig_id)
        if os.path.exists('ligands/unique/{}_lig.mae'.format(lig_id)): continue
        if os.path.exists('ligands/duplicate/{}_lig.mae'.format(lig_id)): continue
        if not os.path.exists(lig_path): continue

        if len(duplicates.keys()) == 0:
            for u_lig in [l.split('_')[0] for l in os.listdir('ligands/unique')]:
                duplicates[u_lig] = [u_lig]
                lig_st[u_lig] = StructureReader('ligands/unique/{}_lig.mae'.format(u_lig)).next()

        new_ligs_seen += 1
        #try:
        st_list = [st for st in StructureReader(lig_path)] # all generated epik states
        charges = [st.formal_charge for st in st_list]

        if len(st_list) == 0:
            os.system('rm -rf ligands/prepared_ligands/{}_lig'.format(lig_id))
            print 'check lig processing', lig_path
            continue
        if len(st_list) == 1 or lig_id not in ligs: 
            lig_st[lig_id] = st_list[0]
        elif ligs[lig_id].target_prot in ['P14416', 'P61169'] and charges[0] == 0 and charges[1] == 1:
            lig_st[lig_id] = st_list[1]#StructureReader(lig_path).next()
        else:
            lig_st[lig_id] = st_list[0]        

        # put all drugs in the unique folder unconditionally
        if lig_id in drugs:
            print 'drug found',lig_id,l
            st_wr = StructureWriter('ligands/unique/{}_lig.mae'.format(lig_id))
            st_wr.append(lig_st[lig_id])
            st_wr.close()
            new_ligs_seen -= 1
            continue
        #print 'hmm', lig_id, drugs.keys()[:5]
        #return        
        for lig_id2 in duplicates:
            if lig_st[lig_id].isEquivalent(lig_st[lig_id2], True):
                duplicates[lig_id2].append(lig_id)
                break
        else:
            duplicates[lig_id] = [lig_id]

    if new_ligs_seen > 0: print 'found {} new ligands'.format(new_ligs_seen)
    
    for lig_id, dups in duplicates.items():
        pdb_ligs = sorted([l for l in dups if l[:6] != 'CHEMBL'])
        specified_stereo = [l for l in dups if l in ligs and ligs[l].valid_stereo]
        
        ki_sort = lambda x: ligs[x].ki

        if len(pdb_ligs) > 0:
            keep = pdb_ligs[0] # alphabetically first non-chembl ligand
        elif len(specified_stereo) > 0:
            specified_stereo.sort(key=ki_sort)
            keep = specified_stereo[0] # best ki ligand with specified stereochemistry
        else:
            dups.sort(key=ki_sort)
            keep = dups[0] # best ki ligand
        
        for l in dups:
            path_write = 'ligands/duplicate/{}_lig.mae'.format(l)
            path_delete = 'ligands/unique/{}_lig.mae'.format(l)
            if l == keep:
                path_write, path_delete = path_delete, path_write
            os.system('rm -f {}'.format(path_delete))
            if os.path.exists(path_write): continue
            st_wr = StructureWriter(path_write)
            st_wr.append(lig_st[l])
            st_wr.close()




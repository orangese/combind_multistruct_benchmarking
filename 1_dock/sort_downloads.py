import os
import sys

from datetime import date

from schrodinger.structure import StructureReader, StructureWriter
from schrodinger.structutils.analyze import evaluate_asl

class PDB:
    def __init__(self, pdb, chain, title, date, lig, prot):
        self.pdb = pdb
        self.chain = chain
        self.title = title
        self.date = date
        self.lig = lig
        self.prot = prot

# taken from schrodinger "solvent" asl definition
solvent = ['hoh', 'spc', 't4p', 't3p', 't5p', 't4pe', 'dod', 'GOL', 'EDO', 'DMS', 'PEG', 'MPD', 
           'PG4', 'BME', 'PGE', 'IPA', '1PE', 'EOH', 'P6G', 'DMF', 'MOH', 'POL', '2PE', 'PE4',
           'ETA', '15P', '12P', 'P33', 'PE5', 'TBU', 'CCN', 'CXE', '7PE', 'NME', 'PE3', 'PE6', 'NHE',
           'ACN', 'P4C', 'PE8', 'NEH', 'XPE', '1BO', 'N8E', 'DMN', 'CE1', 'PYE', 'C10', 'SBT', 'MB3']
solvent = [x.lower() for x in solvent]

ignore_titles = ['mutant', 'mutation', 't877a', 'w741l', 'cryptic']#, 'fragment']#, 'allosteric']
orthosteric = ['dht','tes'] # ar
            
allosteric = ['8vb','8z5','z24','imw', # plk1
              '2an', # cdk2
              'ld2'] # mr
ignore_ligands = solvent + allosteric + ['', 'dtt', 'epe', 'mes', 'tla', 'tar', 'nmm', 
                 'eu','au','edt', 'ola','olb','olc','css', 'cso', 'clr', 'pg0', 'srt', 'cxs', 'tpo', 'sog',
                 'ccs','cme','ben','glc','mli','ocs','aly','sgm','ptr','csd','kcx','ste','cit','iod','hto',
                 'hez','jzr','cps','bog','acp','anp','adp','atp','ags']

ignore_pdb = ['1h00','1h01','1h07','1h08', # cdk2 2 small molecules on top of each other?
              '5tlx','5kcf','5kct','5kra','5tlp', # era multiple small molecules per structure
              '4h71'] # plk1 allosteric site

ignore_methods = ['solution nmr']

def read_information(log):
    pdb_csv = [f for f in os.listdir('structures/downloads') 
                   if f.split('.')[-1] == 'csv' and f[0] != '.']

    # load in csv file
    structs = {}
    for csv in pdb_csv:
        log.write("Begin processing CSV file {}.\n".format(csv))
        uniprot = csv.split('.')[0].lower()

        with open('structures/downloads/{}'.format(csv)) as f:
            for line in f:
                line = line.lower()
                if len(line) <= 1: continue
                    
                # Read header
                if line[:6] == 'pdb id':
                    line = [s.strip() for s in line.split(',')]
                    pdb_i = 0
                    chain_i = line.index('chain id')
                    title_i = line.index('structure title')
                    date_i = line.index('rel. date')
                    lig_i = line.index('ligand id')
                    ligmw_i = line.index('ligand mw')
                    prot_i = line.index('uniprot acc')
                    method_i = line.index('exp. method')
                    continue
                line = [s.strip('"') for s in line.split('","')]
                u_list = [p.strip('"\n') for p in [x.strip() for x in line[prot_i].strip().split(',')]]                 
                y,m,d = line[date_i].split('-')
                st_date = date(int(y),int(m),int(d))
                skipping = "Skipping PDB {}, chain {}".format(line[pdb_i], line[chain_i])

                if line[pdb_i] not in structs: structs[line[pdb_i]] = []
                
                if line[pdb_i] in ignore_pdb:
                    log.write("{}: PDB ID {} in ignore_pdb.\n".format(skpping, line[pdb_i]))
                    continue
                if line[lig_i] in ignore_ligands:
                    log.write("{}: Ligand {} in ignore_ligands.\n".format(skipping, line[lig_i]))
                    continue
                if line[method_i] in ignore_methods:
                    log.write("{}: Method {} in n ignore_ligands.\n".format(skipping, line[method_i]))
                    continue
                if uniprot not in u_list:
                    log.write("{}: Uniprots {} not target uniprot.\n".format(skipping, ','.join(u_list)))
                    continue
                if any(x in line[title_i] for x in ignore_titles):
                    log.write("{}: Titles {} in ignore_titles.\n".format(skipping, line[title_i]))
                    continue
                if float(line[ligmw_i]) < 100 or float(line[ligmw_i]) > 1000:
                    log.write("{}: Ligand of MW {} not right size.\n".format(skipping, line[ligmw_i]))
                    continue
                
                structs[line[pdb_i]] += [PDB(line[pdb_i], line[chain_i], 
                                             line[title_i], st_date, line[lig_i], uniprot)]
        log.write("End processing CSV file {}.\n".format(csv))
    return structs
    
def get_ligs(structs, log):
    log.write('Begin dismabiguating ligands.\n')
    ligs = {}
    for pdbid, pdblist in sorted(structs.items(),key=lambda x:x[0]):

        if len(set([p.lig for p in pdblist])) != 1:
            log.write("Ignoring PDB {}: Unable to choose ligand from {}.\n".format(pdbid, set([p.lig for p in pdblist])))
            continue
        if pdblist[0].lig not in ligs: ligs[pdblist[0].lig] = []    
        first_chain = sorted(pdblist, key=lambda x:x.chain)[0] 
        ligs[pdblist[0].lig].append(first_chain)
    log.write('End disambiguating ligands')
    return ligs

def write_files(ligs, log):
    log.write('Begin file writing.\n')
    for lname, l_copies in sorted(ligs.items(),key=lambda x:x[1][0].pdb):
        if any(os.path.exists('structures/raw_files/{}_lig.mae'.format(l.pdb.upper())) for l in l_copies):
            log.write("Ignoring PDBs {}: Ligand already present in raw_files\n".format(','.join(l_copies)))

        l_copies.sort(key=lambda x: x.date)

        # find the first non-covalent copy
        for l in l_copies: 
            prot_st = StructureReader('structures/downloads/{}.pdb'.format(l.pdb)).next()
            lig_st = None
            for res in prot_st.residue:
                if res.pdbres.strip().lower() == lname and res.chain.strip().lower() == l.chain:
                    molnum = res.molecule_number
                    if len([m.residue for m in prot_st.molecule if m.number == molnum][0]) == 1:
                        lig_st = prot_st.extract([a.index for a in res.atom])
                    break
            if lig_st is None:
                log.write("Ignoring PDB {}: Ligand covalently bound\n".format(l.pdb))
                continue
            
            for pdb in l_copies[l_copies.index(l)+1:]:
                log.write("Ignoring PDB {}: Ligand redundant with {}\n".format(pdb.pdb, l.pdb))

            log.write("Writing PDB {}.\n".format(l.pdb))
            to_delete = evaluate_asl(prot_st, 'solvent or res.pt {}'.format(lname.upper()))
            prot_st.deleteAtoms(to_delete)
                
            prot_wr = StructureWriter('structures/raw_files/{}_prot.mae'.format(l.pdb.upper()))
            prot_wr.append(prot_st)
            prot_wr.close()
            
            lig_wr = StructureWriter('structures/raw_files/{}_lig.mae'.format(l.pdb.upper()))
            lig_wr.append(lig_st)
            lig_wr.close()
            break
    log.write('End file writing.\n')


def sort_downloads():
    assert not os.path.exists('pdb_processing.log'), 'PDB processing has already been run'
    with open('pdb_processing.log', 'a') as log:
        log.write('Manual processing\n')
        log.write('Begin automated processing logging\n')
        
        structs = read_information(log)
        log.write("Found at least 1 valid entry for {} of {} PDBs.\n".format(sum(len(v) for v in structs.values()), len(structs)))

        ligs = get_ligs(structs, log)

        os.system('mkdir -p structures/raw_files')
        write_files(ligs, log)
        log.write('End automated processing logging\n')
    return ligs

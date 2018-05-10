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

ignore_title = ['mutant', 'mutation', 't877a', 'w741l', 'cryptic']#, 'fragment']#, 'allosteric']
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

def sort_downloads():
    
    pdb_csv = [f for f in os.listdir('structures/downloads') 
               if f.split('.')[-1] == 'csv' and f[0] != '.']

    # load in csv file
    for csv in pdb_csv:
        uniprot = csv.split('.')[0].lower()
        structs = {}
        with open('structures/downloads/{}'.format(csv)) as f:
            for line in f:
                line = line.lower()# [s.strip('"').lower() for s in line.strip().split(',')]

                if len(line) <= 1: continue
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
                 
                if line[pdb_i] in ignore_pdb: continue
                if line[lig_i] in ignore_ligands: continue
                #print line[pdb_i],line[lig_i]
                if line[method_i] == 'solution nmr': continue
                u_list = [p.strip('"\n') for p in [x.strip() for x in line[prot_i].strip().split(',')]]
                
                if uniprot not in u_list: continue
                if True in [x in line[title_i] for x in ignore_title]: continue
                if float(line[ligmw_i]) < 100 or float(line[ligmw_i]) > 1000: continue
                
                y,m,d = line[date_i].split('-')
                st_date = date(int(y),int(m),int(d))
                
                if line[pdb_i] not in structs: structs[line[pdb_i]] = []
                structs[line[pdb_i]].append(PDB(line[pdb_i], line[chain_i], 
                    line[title_i], st_date, line[lig_i], uniprot))

    # sort by ligand
    ligs = {}
    for pdbid, pdblist in sorted(structs.items(),key=lambda x:x[0]):
        lig_name = None
        if len(set([p.lig for p in pdblist])) == 1: 
            lig_name = pdblist[0].lig
        
        if lig_name is not None:
            first_chain = sorted(pdblist, key=lambda x:x.chain)[0] 
            if lig_name not in ligs:# or ligs[lig_name].date > first_chain.date:
                ligs[lig_name] = []#first_chain
            ligs[lig_name].append(first_chain)

        else:
            for chain in [p.chain for p in pdblist]:
                for lig in [p.lig for p in pdblist if p.chain == chain]:
                    if lig in orthosteric:
                        break
                else:
                    print pdbid, [(p.chain, p.lig) for p in pdblist]

    # write out ligand/protein structures
    #print len(ligs), 'ligands found'
    os.system('mkdir -p structures/raw_files')
    for lname, l_copies in sorted(ligs.items(),key=lambda x:x[1][0].pdb):
        if True not in [os.path.exists('structures/raw_files/{}_lig.mae'.format(l.pdb.upper())) for l in l_copies]:
            #print lname, l.pdb, l.date
            #print l.title

            l_copies.sort(key=lambda x: x.date)

            for l in l_copies: # find the first non-covalent copy
                #print l.lig, l.pdb, l.date
                prot_st = StructureReader('structures/downloads/{}.pdb'.format(l.pdb)).next()
                lig_st = None
                for res in prot_st.residue:
                    if res.pdbres.strip().lower() == lname and res.chain.strip().lower() == l.chain:
                        molnum = res.molecule_number
                        if len([m.residue for m in prot_st.molecule if m.number == molnum][0]) == 1:
                            lig_st = prot_st.extract([a.index for a in res.atom])
                        else:
                            pass#print lname, l.pdb, 'covalent ligand'
                        break
                if lig_st is None:
                    #print lname, 'not found'
                    continue
                #break#continue
                to_delete = evaluate_asl(prot_st, 'solvent or res.pt {}'.format(lname.upper()))
                prot_st.deleteAtoms(to_delete)
            
                prot_wr = StructureWriter('structures/raw_files/{}_prot.mae'.format(l.pdb.upper()))
                prot_wr.append(prot_st)
                prot_wr.close()
            
                lig_wr = StructureWriter('structures/raw_files/{}_lig.mae'.format(l.pdb.upper()))
                lig_wr.append(lig_st)
                lig_wr.close()
                break            

    return ligs









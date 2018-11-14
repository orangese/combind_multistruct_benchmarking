import os
import sys
sys.path.append('../../1_dock')
sys.path.append('../../2_ifp')
sys.path.append('../../3_analyze')
from containers import Dataset
from shared_paths import shared_paths

template = '/scratch/PI/rondror/combind/bpp_data/{0:}/docking/mininplace/{1:}-to-{2:}/{1:}-to-{2:}.rept'
shared_paths['docking'] = 'XP'
fname = '/scratch/PI/rondror/combind/bpp_outputs/glide_performance_{}.tsv'.format(shared_paths['docking'])
with open(fname, 'w') as out:
    for protein in os.listdir(shared_paths['data']):
        print(protein)
        if protein[0] == '.': continue
        data = Dataset(shared_paths, [protein])
        lm = data.proteins[protein].lm
        print(lm.pdb)

        data.load({protein: lm.pdb}, {protein: [lm.st]}, load_mcss = False, load_fp = False)

        docking =  data.proteins[protein].docking[lm.st]
        for lig_name, ligand in docking.ligands.items():
            mininplace = template.format(protein, lig_name, lm.st)
            if os.path.exists('/'.join(mininplace.split('/')[:-1])):
                gscore, emodel, rmsd = 'inf', 'inf', 'inf'
            else:
                gscore, emodel, rmsd = '.', '.', '.'
            if os.path.exists(mininplace):
                with open(mininplace) as fp:
                    if 'NO ACCEPTABLE' not in fp.readline():
                        for line in fp:
                            if line[0] == '=': break
                        line = fp.readline()
                        rank, lig, score, gscore = line.strip().split()[:4]
                        emodel = line.split()[13]
                        assert rank == '1', (rank, protein, lig_name)
                        assert lig == lig_name
                        try:
                            with open('/'.join(mininplace.split('/')[:-1])+'/rmsd.csv') as rmsd_file:
                                rmsd_file.readline()
                                rmsd = rmsd_file.readline().split(',')[3].strip('"')
                        except:
                            pass

            out.write('\t'.join([
                protein,
                lm.st,
                lig_name,
                ','.join(map(lambda x: str(x.rmsd), ligand.poses)),
                ','.join(map(lambda x: str(x.gscore), ligand.poses)),
                ','.join(map(lambda x: str(x.emodel), ligand.poses)),
                gscore,
                emodel,
                rmsd
            ]) + '\n')

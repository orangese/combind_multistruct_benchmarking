import os
import sys

data = '/scratch/PI/rondror/jbelk'

dock_st = {'B1AR':'2VT4','AR_final':'3B67','TRPV1':'3J5Q'}
all_ligs = {p:[l.split('.')[0] for l in os.listdir('{}/{}/unique_ligands'.format(data, p))
               if l != st+'_lig.mae'] for p,st in dock_st.items()}

#out_dir = 'scores3'
#all_lam = [i/10.0 for i in range(1,11)]
#all_w10 = [i/100.0 for i in range(6)]
#all_n = [i for i in range(21)]
#all_w56 = [(0,0),(0,1),(1,0)]
# ---------------------------------
out_dir = 'scores5'
all_lam = [i/100.0 for i in range(5,16)]
all_w10 = [i/200.0 for i in range(4)]
all_n = [i for i in range(51)]
all_w56 = [(0,0),(0,1),(1,0)]

#print all_lam
#print all_w10
#exit()

for prot in sorted(dock_st.keys()):
    #if prot == 'AR_final': continue
    for lig in sorted(all_ligs[prot]):
        print prot, lig
        #if lig != '5IS0_lig': continue
        for x in enumerate(all_lam):
            for y in enumerate(all_w10):
                if y[0] != 0: continue
                #for z in enumerate(all_n):
                for m in enumerate(all_w56):
                    if m[0] != 0: continue
                        #all_settings.append((query, x, y, z, m))
                    for z in all_n:
                        f_name = '{}-{}-{}-{}.txt'.format(x[0], y[0], z, m[0])
                        f_path = '{}/{}/{}/{}/{}'.format(data, prot, out_dir, lig, f_name)
                        if not os.path.exists(f_path):
                        
                        #if os.path.exists(f_path): continue
                        #print prot, lig, x[0], y[0], z[0], m[0]
                            os.system('sbatch -p owners -t 2:30:00 --job-name={}-{}-{}-{} score.sh {} {} {} {} {}'.format(lig, x[0], y[0], m[0], prot, lig, x[0], y[0], m[0]))
                            break
                        else: continue#print 'done', lig, x, y, m
        #if lig not in submit: continue
        #print prot, lig
        #os.system('sbatch --tasks=6 -p owners --job-name={} score.sh {} {}'.format(lig, prot, lig))
        #break
    #break

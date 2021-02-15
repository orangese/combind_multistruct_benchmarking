import pandas as pd
import os
import subprocess
index = '/oak/stanford/groups/rondror/projects/ligand-docking/pdbbind_2019/refined-set-raw/index/INDEX_refined_name.2019'
root = '/oak/stanford/groups/rondror/projects/ligand-docking/pdbbind_2019/chembl'

uniprots = set()
with open(index) as fp:
	for line in fp:
		if line[0] == '#': continue
		uniprots.add(line.split()[2])

print('Read {} uniprots.'.format(len(uniprots)))

for uniprot in uniprots:
	fname = '{}/{}_all.csv'.format(root, uniprot)
	log = '{}/{}_all.log'.format(root, uniprot)
	if not os.path.exists(fname):
		cmd = ('python chembl/chembl.py --homologous --ambiguous-stereo '
		       '--output-fname {} {}')
		cmd = cmd.format(fname, uniprot)

		sbatch = 'sbatch -p owners -t 04:00:00 -J chembl --output={} --wrap="{}"'
		sbatch = sbatch.format(log, cmd)

		subprocess.run(sbatch, shell=True)

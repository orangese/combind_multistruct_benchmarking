from pymol import cmd


def helpers(fname, query):
	prot = fname.split('/')[0]

	with open(fname) as fp:
		for line in fp:
			q, h = line.strip().split(':')
			if q == query:
				h = h.split(',')

	for ligand in [query] + h:
		cmd.load('{0:}/ligands/prepared_ligands/{1:}/{1:}.mae'.format(prot, ligand))

	cmd.alignto(query)


cmd.extend('helpers', helpers)
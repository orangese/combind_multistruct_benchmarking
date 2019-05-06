import sys
import os
from glob import glob
from schrodinger.structure import StructureReader, SmilesWriter
from multiprocessing.dummy import Pool as ThreadPool 

def setup():
	"""
	Creates directory ligands/from_smiles
	and fills it with directories for each pdb ligand
	containing the raw ligand structure
	"""
	for protein in glob('*'):
		os.system('mkdir -p {}/ligands/from_smiles')
		for ligand in sorted(glob('{}/structures/raw_files/*_lig.mae'.format(protein)))[:20]:
			name = ligand.split('/')[-1].split('_')[0]
			print(name)
			d = '{}/ligands/from_smiles/{}'.format(protein, name)
			os.system('mkdir -p {}'.format(d))
			os.system('cp {} {}'.format(ligand, d))

def prep(path):
	print(path)
	os.chdir(path)
	pdb = path.split('/')[-1]
	if os.path.exists('{}_lig_prep.mae'.format(pdb)): return
	os.system('{0}/utilities/prepwizard -WAIT -noepik -rehtreat -noprotassign -noimpref '
		      '{1}_lig.mae {1}_lig_prep.mae'.format(os.environ['SCHRODINGER'], pdb))


def pdb_to_smiles(ligand, smiles_fname):
	"""
	Converts processed ligands in MAE format into SMILES strings.
	
	Inputs:
		ligands ([str, ...])
		smiles_fname (str)

	Returns None
	"""
	if os.path.exists(ligand):
		return
	

	with SmilesWriter(smiles_fname) as out:
		with StructureReader(ligand) as sts:
			for st in sts:
				out.append(st)

def all_pdb_to_smiles():
	mae =  list(glob('*/ligands/from_smiles/*/*_lig_prep.mae'))
	smi = [m.replace('_lig_prep.mae', '.smi') for m in mae]

	pool = ThreadPool(1)
	pool.map(pdb_to_smiles, zip(mae, smi))

def ligprep():
	for path in glob('/scratch/PI/rondror/combind/bpp_data/*/ligands/from_smiles/*'):
		print(path)
		os.chdir(path)
		pdb = path.split('/')[-1]
		os.system('{0}/utilities/ligprep -WAIT -epik '
		          '{1}.smi {1}.mae'.format(os.environ['SCHRODINGER'], pdb))


if sys.argv[1] == 'setup':
	setup()
elif sys.argv[1] == 'prep':
	paths = glob('/scratch/PI/rondror/combind/bpp_data/*/ligands/from_smiles/*')
	pool = ThreadPool(1)
	pool.map(prep, paths)
elif sys.argv[1] == 'smi':
	all_pdb_to_smiles(*sys.argv[2:])
"""
Wrapper to use Schrodinger shape_screen for virtual screening.

Given
	- Template ligand pose, generally an XTAL conformation "template.maegz"
	- Known binders "known.maegz"
	- Screening library "ligands.maegz"

# Write known binders to one file.
python shape_screen.py extract --best paths/to/known/*.maegz known.maegz

# Align known molecules to the template molecule, often an XTAL conformation.
python shape_screen.py screen template.maegz known.maegz

# Merge aligned binders into a single entry.
python shape_screen.py merge known_align.maegz known_align_merged.maegz

# Screen library against the aligned, merged known binders
python shape_screen.py screen known_align_merged.maegz ligands.maegz
"""

import os
import click
from schrodinger.structure import StructureReader, StructureWriter

@click.group()
def main():
	pass

@main.command()
@click.argument('merged')
@click.argument('paths', nargs=-1)
@click.option('--best', is_flag=True)
def extract(merged, paths, best):
	"""
	merged (mae file path): where to write combined set of ligands
	paths [mae file path, ...]: paths to mae files for ligands.
	"""
	with StructureWriter(merged) as writer:
		for path in paths:
			with StructureReader(path) as reader:
				for st in reader:
					writer.append(st)
					if best: break

@main.command()
@click.argument('input_mae')
@click.argument('output_mae')
def merge(input_mae, output_mae):
	"""
	merged (mae file path): where to write combined set of ligands
	paths [mae file path, ...]: paths to mae files for ligands.
	"""
	with StructureReader(input_mae) as sts:
		st = next(sts)
		for _st in sts:
			st = st.merge(_st)
		st.write(output_mae)

@main.command()
@click.argument('template')
@click.argument('ligands')
@click.option('--aligned', default='{ligands}-to-{template}')
def screen(template, ligands, aligned):
	aligned = aligned.format(template=template.split('/')[-1].split('.')[0],
		                     ligands=ligands.split('/')[-1].split('.')[0])

	out_mae = aligned + '_align.maegz'
	out_csv = aligned + '_align.csv'

	if not os.path.exists(out_mae):
		cmd = '$SCHRODINGER/shape_screen -shape {template} -screen {ligands} -WAIT -flex -pharm -JOB {aligned}'
		cmd = cmd.format(template=template, ligands=ligands, aligned=aligned)
		print(cmd)
		os.system(cmd)

	with StructureReader(out_mae) as sts, open(out_csv, 'w') as fp:
		fp.write('ID,variant,score\n')
		for st in sts:
			fp.write('{},{},{}\n'.format(st.title, st.property['s_lp_Variant'], st.property['r_phase_Shape_Sim']))

main()

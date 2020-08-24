"""
Wrapper to use Schrodinger shape_screen for virtual screening.

Given
	- Template ligand pose, generally an XTAL conformation
	- Known binders
	- Screening library
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
def merge(merged, paths, best):
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
@click.argument('template')
@click.argument('ligands')
@click.option('--aligned', default='')
def screen(template, ligands, aligned):
	if not aligned:
		aligned = ligands.split('/')[-1].split('.')[0]
		aligned += '-to-'
		aligned += template.split('/')[-1].split('.')[0]

	cmd = '$SCHRODINGER/shape_screen -shape {template} -screen {ligands} -WAIT -flex -pharm -JOB {aligned}'
	cmd = cmd.format(template=template, ligands=ligands, aligned=aligned)
	print(cmd)
	os.system(cmd)

main()

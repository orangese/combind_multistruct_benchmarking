import sys
import pandas as pd
from schrodinger import structure  
import os
import glob
import click

def get_glide_file_for_ligand(ligand, file_names):
    for file_name in file_names:
        if ligand in file_name:
            return file_name
    return -1

@click.command()
@click.argument('root_directory', type=click.Path(exists=True))
@click.argument('out_directory', type=click.Path(exists=True))
@click.argument('system')
@click.option('--docking-version', default='confgen_es4')
@click.option('--scores', default='scores/pdb.sc')
def main(root_directory, out_directory, system, docking_version, scores):
	struct_file_paths = glob.glob("{}/{}/docking/{}/*/*maegz".format(root_directory, system, docking_version))
	df_glide_scores = pd.read_csv("{}/{}/{}".format(root_directory, system, scores))

	for idx, row in df_glide_scores.iterrows():
	    glide_file = get_glide_file_for_ligand(row["lig"], struct_file_paths)
	    if glide_file == -1:
	        print ('missing file for '+str(row['lig']))
	        continue
	    
	    with structure.StructureReader(glide_file) as reader:
	        #indexing from -1 so that first ligand pose has index = 0!
	        i = -1
	        for pose in reader:
	            if i == -1:
	                pose.write("{}/protein.mae".format(out_directory))
	            
	            if i == int(row["combind_rank"]):
	            	pose.write("{}/{}_top_combind.mae".format(out_directory, row["lig"])) 
	            
	            if i == 0:
	                pose.write("{}/{}top_ten_glide.mae".format(out_directory, row['lig']))
	            if 0 < i < 10:
	            	pose.append("{}/{}top_ten_glide.mae".format(out_directory, row['lig']))
	            
	            if i > 10 and i >= int(row["combind_rank"]):
	                break
	            i += 1

main()
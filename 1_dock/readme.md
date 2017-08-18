Make a new dataset, process all files, and get docking results.

## 1) Get structures from the Protein Data Bank (PDB)

Go to `www.rcsb.org` and search for a receptor. Download all files (button at the top of the results page). The Download Tool should open. In the first section, Dowload Coordinates and Experimental Data, select two options: PDB Coordinates and Uncompressed. Deselect all other options. Launch the download.

Move these files into `/scratch/PI/rondror/docking_data/<receptor name>/raw_pdbs/`. You may need to create new directories:

`cd /scratch/PI/rondror/docking_data/`

`mkdir <receptor name>`

`mkdir <receptor name>/raw_pdbs`

## 2) Process files and dock ligands

Remove alternate conformations from the pdb files:

`sbatch xdock.sh -r <receptor>`

Strip waters and align structures:

`sbatch xdock.sh -s <receptor>`

Add hydrogens and missing sidechains:

`sbatch xdock.sh -p <receptor>`

Extract ligands:

`sbatch xdock.sh -l <receptor>`

Define the receptor region the ligand is docked into (the Glide grid):

`sbatch xdock.sh -g <receptor>`

Dock all the ligands into all the grids:

`sbatch xdock.sh -x <receptor>`

For example:

`sbatch xdock.sh -x B1AR`

## 3) Evaluate docking quality

Before proceeding to fingerprinting, make sure to check the docking results using the heatmap plot in `notebooks`. If there is not at least one receptor that generates good poses (at least one pose with an rmsd < 2A) for most ligands, do not continue! Either reconfigure glide to get better docking results or choose another receptor.

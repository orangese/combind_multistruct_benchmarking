Make a new dataset, process all files, and get docking results.

## 1) Get structures from the Protein Data Bank (PDB)

Go to `www.rcsb.org` and search for a receptor. Download all files (button at the top of the results page). The Download Tool should open. In the first section, Dowload Coordinates and Experimental Data, select two options: PDB Coordinates and Uncompressed. Deselect all other options. Launch the download.

Move these files into `/scratch/PI/rondror/docking_data/<receptor name>/raw_pdbs/`. You may need to create new directories:

`cd /scratch/PI/rondror/docking_data/`

`mkdir <receptor name>`

`mkdir <receptor name>/raw_pdbs`

## 2) Process files and dock ligands automatically

The following can be done in one step with the command:

`./main.py -splgd <receptor name>`.

For example:

`./main.py -splgd B2AR`

There are five steps: Strip, Process, Ligands, Grids, and Dock. Each can be run independently with, for example:

`./main.py -s <receptor name>`.

## 3) Compute RMSDs for all poses

`./compute_rmsds.py <receptor name>`

RMSDs will be output to `/scratch/PI/rondror/docking_data/<receptor name>/rmsd.csv`.

## 4) Evaluate docking quality

Before proceeding to fingerprinting, make sure to check the docking results using the heatmap plot in `notebooks`. If there is not at least one receptor that generates good poses (at least one pose with an rmsd < 2A) for most ligands, do not continue! Either reconfigure glide to get better docking results or choose another receptor.

How to make a new dataset, process all files, and get docking results.

## 1) Get structures from the Protein Data Bank (PDB)

Go to `www.rcsb.org` and search for a receptor. Download all files (button at the top of the results page). The Download Tool should open. In the first section, Dowload Coordinates and Experimental Data, select two options: PDB Coordinates and Uncompressed. Deselect all other options. Launch the download.

Move these files into `/scratch/PI/rondror/docking_data/<receptor name>/raw_pdbs/`. You may need to create new directories:

`cd /scratch/PI/rondror/docking_data/`

`mkdir <receptor name>`

`mkdir <receptor name>/raw_pdbs`

## 2) Process files and dock ligands automatically

The following can be done in one step with the command:

`./main.py -splgd <receptor name>`.

There are five steps: Strip, Process, Ligands, Grids, and Dock. Each can be run independently with, for example:

`./main.py -s <receptor name>`.

## 3) Compute RMSDs for all poses



## 4) Evaluate docking quality



This folder contains scripts to automatically process files for docking. There are two necessary inputs: a list of chembl ligands and the structure files.

0. CHEMBL input:
- <protein>/chembl
    - find the correct target on chembl (e.g., D2 is CHEMBL217)
    - download either Ki or IC50 data. I usually pick Ki data unless there is way more IC50 data (e.g., for some kinases there is usually more IC50 data). 
    - it is fine to download ki data for multiple organisms (e.g., rat and human) and put them in this folder
    - the processing step will merge multiple files
    - see D2R/chembl for an example

1. Structure input:
- <protein>/structures/raw_files
    - <pdb>_lig.mae and <pdb>_prot.mae for each <pdb> of interest
run `$SCHRODINGER/run prepare_all.py <prot1> <prot2> ...
(Under development/ignore for now):
Automatic structure input:
- <protein>/structures/downloads
    - do a pdb search by uniprot id
    - download all .pdb files
    - from the pdb search results page, download structure report
        - select "Customizable Table" from drop down menu
        - in the new page "Custom Tabular Report", select: (1) the entire "Structure Summary" section, (2) "Ligand ID" and "Ligand Molecular Weight" under "Ligands", and (3) "Uniprot Accessions" under "Sequence Clusters".
        - the point of this report is to automatically parse relevant ligands out of the structures. This is very helpful when dealing with proteins with many structures.
        - create the report and download the csv. Put the .csv in <protein>/structures/downloads with the structures.
        - you must save this file as <uniprot>.csv and the <uniprot> must match the <uniprot> of the downloaded structures (this is used when sorting out the downloaded files).

1. Prepare structures
- process_structs() merges protein and ligand and prepares them together at pH 7.0
- align_structs() aligns all the structures to a single template structure for use in computing the ligand rmsd
    - align_structs() also attempts to renumber the proteins to be consistent but this does not always work
- sort_files() separates the proteins and ligands back into separate files (based on which molecule was present in which file in raw_files)
- make_grids() makes a grid for each protein. the centroid of the co-crystallized ligand defines the center of the grid.
- dock() docks all pairs of ligands with structures. view this output in a notebook (i.e. notebooks/1_view_glide_results) to pick a structure to proceed with for all the chembl ligands/future steps. Put that structure in the `grids` dictionary.
2. Prepare ligands
- reads in all ligands including chembl ligands and structure ligands
- prepares them at pH 7.0, removes extraneous ions/other things from chembl ligands (`desalt`), and removes duplicate ligands
- init_mcss() computes the mcss of each chembl ligand to each structure ligand for use in pick_helpers()
- chembl ligands with molecular weight > 1000 or Ki > 1 uM are ignored.
3. Decide what ligands to use
- pick helper (chembl) ligands for each query based on desired criteria. outputs these as text files to <prot>/chembl/helpers
- compute the mcss of all pairs of ligands needed for each text file
4. Dock, fingerprint, and mcss all ligands
- to whatever protein structure is specified in `grids`

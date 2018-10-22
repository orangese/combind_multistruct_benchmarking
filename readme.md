# ComBind

ComBind scores small molecule ligand binding poses
by combining data-driven and physics-based methods.

## Organization

The code base is seperated into subdirectories by task.
Each directory contains its own readme highlighting
how the code works, its purpose, known issues, and plans
for future development (Some of these documents
are still under development). A general overview of the
workflow is below.

There are currently no automated tests for this project,
however, as examples where the code does not behave as
expected are identified they should be added to the test
directory in each directory. The hope is that we will only
have to manually discover each failure mode once.

Before the code is run, you should source setup_combind.

Inside the code global constants are stored in shared_paths.py.k
You should never pass the shared_paths dictionary
as an argument, as this is unncessary and confusing.

## Workflow Overview

0. Protein and Ligand Curation and Processing

   __Input:__ 

1. Docking
   
   __Input:__ Ligand/receptor complexes downloaded from the Protein Data Bank and ligands from CHEMBL.
   
   __Output:__ Docking results (a set of poses) for each ligand (cross-)docked into each receptor.

2a. Fingerprinting
   
   __Input:__ Docking results and crystal structures.
   
   __Output:__ An interaction fingerprint for each pose.

2b. Maximum Common Substructure (mcss) Computation

   __Input:__ Docking results.
   
   __Output:__ RMSDs between MCSSs for all ligand pairs that share a common scaffold.

3. Scoring Function
   
   __Input:__ Interaction fingerprints for poses for several ligands.
   
   __Output:__ Score for set of poses.

4. Optimization

   __Input:__ Interaction fingerprints for all poses of several ligands. ComBind scoring function.
   
   __Output:__ Set of maximum likelihood poses for given ligands.

5. Notebooks
   
   __Input:__ Scores, fingerprints, and rmsds for all poses for all ligands.
   
   __Output:__ Many tools to visualize our performance and track our progress.

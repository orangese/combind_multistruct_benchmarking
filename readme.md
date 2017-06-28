## Overview

The goal of combind is data-driven pose prediction for small molecules.

## Project Structure

Project status and history can be found in the issues, projects, pull requests, and wiki of this repository.

For information about previous versions of combind, please see the readme in `drorlab/sar_workflow`.

## Running the Code

Our code is organized into four steps:

1. Docking
   
   __Input:__ Ligand/receptor complexes downloaded from the Protein Data Bank.
   
   __Output:__ Docking results (a set of poses) for each ligand (cross-)docked into each receptor.

2. Fingerprinting
   
   __Input:__ Docking results and crystal structures.
   
   __Output:__ An interaction fingerprint for each pose.

3. Scoring
   
   __Input:__ Interaction fingerprints for all poses for all ligands.
   
   __Output:__ Ranked poses for each ligand.

4. Notebooks
   
   __Input:__ Scores, fingerprints, and rmsds for all poses for all ligands.
   
   __Output:__ Many tools to visualize our performance and track our progress.

See the readme in each directory for the commands to run to execute each step.

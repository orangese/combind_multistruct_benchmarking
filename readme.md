# ComBind

ComBind scores small molecule ligand binding poses
by combining data-driven and physics-based methods.

## Creating a Schrodinger virtual environment

ComBind makes use of the schrodinger python api, which is only accessible by
using their interpretter. Creating a virtual environment that makes their
interpretter the default python interpretter is the simplest way to do this.
By following the below directions you'll be able to use their python api
both in scripts and in jupyter notebooks.

Set the SCHRODINGER environment variable using, for local installation,
`export SCHRODINGER=/opt/schrodinger/suites2019-4` or, for sherlock,
`ml chemistry schrodinger`.


Then create the environment and upgrade the relevant packages.
```
$SCHRODINGER/run schrodinger_virtualenv.py schrodinger.ve
source schrodinger.ve/bin/activate
pip install --upgrade jupyter matplotlib numpy sklearn scipy pandas
```

Run `source schrodinger.ve/bin/activate` to activate the
environment in the future.

## Organization

The code base is seperated into subdirectories by task.
Each directory contains its own readme highlighting
how the code works, its purpose, known issues, and plans
for future development (Some of these documents
are still under development). A general overview of the
workflow is below.

Before the code is run, you should source setup_combind.

## Testing

Tests are located in each directory. They largely fall into two categories
(1) tests of logic of python code and (2) tests of behaviour of schrodinger
code. The former can be run by executing "pytest" from $COMBINDHOME or by
"pytest /path/to/test.py" to run an individual test suite. The latter need
to be run individually and largely do not have automated checking i.e. you
need to manually inspect the output to see if it is what you expect.

Any time you notice that a ligand was not given the right protonation state,
a fingerprint is wrong, proteins are not properly aligned, etc., you should add it
as a test case.

## Workflow Overview

0. Protein and Ligand Curation and Processing

1. Docking
   
   __Input:__ Ligand/receptor complexes downloaded from the Protein Data Bank and ligands from CHEMBL.
   
   __Output:__ Docking results (a set of poses) for each ligand (cross-)docked into each receptor.

2. (a) Fingerprinting
   
   __Input:__ Docking results and crystal structures.
   
   __Output:__ An interaction fingerprint for each pose.

2. (b) Maximum Common Substructure (mcss) Computation

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


## Evaluation Strategy Overview

- __Pose Generation__
  * Fraction of ligands where top pose is correct?
  * Fraction of ligands with any correct pose?
  - __Protein Preparation__
    * Are proteins properly aligned?
    * Do proteins have reasonable protonation states?
  - __Ligand Preparation__
    * Do ligands get reasonable protonation states?
  - __Docking__

- __ComBind Rescoring__
  * Performance improvement over Glide
  - __Scoring Function__
    * How well does the overall scoring function separate native and decoy poses?
    * ... pose pairs? (More tractable to quickly evaluate)
    - __Fingerprinting__
      * Are all observed interactions correctly detected?
    - __Fingerprint Similarity Metric__
      * Are energies similar when merged across different normalization factors
        as compared to when computed for a single value of the normalization factor?
    - __MCSS__
      * How many MCSS do we get?
      * Do substructures tend to overlap in correct poses?
      * Do cases where they don't seem reasonable?
      * How well are native and decoy poses seperated?
    - __Statistics Computation__
      * Are the results smooth?
      * But not too biased?
      * Monotonic?
  - __Optimization__
    * Do we reach the global optima?

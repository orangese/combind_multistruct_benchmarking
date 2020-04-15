# ComBind

ComBind scores small molecule ligand binding poses by combining data-driven
modeling and physics-based docking.

## Installation

Running ComBind requires access to Glide and the Schrodinger Python API.

### Creating a Schrodinger virtual environment

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

## Overview

Running ComBind can be broken into several components: data curation,
data preparation (including docking), featurization of docked poses,
the ComBind scoring itself, and inspection of results.

### Curation of raw data

To produce poses for a particular protein, you'll need to provide at least one
3D structure of the target protein and chemical structures of ligands to dock.

### Data preparation and docking

### Featurization: interaction fingerprints, maximum common substructures

### ComBind Scoring

### Visualization of the results

### Model Training (Optional)

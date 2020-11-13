# ComBind

ComBind scores small molecule ligand binding poses by combining data-driven
modeling and physics-based docking.

Specifically, given the chemical structures of several ligands that can bind
a given target protein, ComBind solves for a set of poses, one per ligand, that
are both highly scored by physics-based docking and display similar interactions
with the target protein. ComBind quantifies this vague notion of "similar" by
considering a diverse training set of protein complexes and computing the
overlap between protein–ligand interactions formed by distinct ligands when
they are in their correct poses, as compared to when they are in randomly
selected poses.

## Predicting poses for known binders

- (ligands.smi) a set of ligands known to bind the target of interest
- (grid.zip) docking grid
- (XTAL_lig.mae) crystal conformation of select ligands, optional
- (stats/*.de) statistics to use in combind scoring
- parameters for pose prediction

Running ComBind can be broken into several components: data curation,
data preparation (including docking), featurization of docked poses,
the ComBind scoring itself, and inspection of results.

See instructuctions for software installation at the bottom of this page.

### Curation of raw data

To produce poses for a particular protein, you'll need to provide at least one
3D structure of the target protein and chemical structures of ligands to dock.

These raw inputs need to be properly stored so that the rest of the pipeline
can recognize them.

The structure(s) should be stored in a directory `structures/raw`.
Each structure should be split into two files `NAME_prot.mae`
and `NAME_lig.mae` containing only the protein and only the reference ligand,
respectively.

Alternatively, if you'd prefer to prepare your structures yourself, save your
prepared files to `structures/proteins` and `structures/ligands`. Moreover,
you could even just begin with a Glide docking grid which you prepared yourself
by placing it in `docking/grids`.

Ligands should be specified in a smiles file with a header line containing at
least the entries "ID" and "SMILES", specifying the ligand name and the ligand
chemical structure.

### Data preparation and docking

```
combind structprep
combind ligprep ligands.smi ligands
combind dock structures/grids/PDBID/PDBID.zip docking ligands/*/*.maegz
```

### Binding pose featurization

```
combind featurize . docking/*/*_pv.maegz
```

### Pose prediction with ComBind
```
combind pose-prediction . poses.csv
```

### Visualizing the results

TODO

## ComBind virtual screening

- (library.smi) library of compounds to screen
- (binders_pv.maegz) poses of known binders to use in combind scoring
- (grid.zip) docking grid
- (stats/*.de) statistics to use in combind scoring
- parameters to use in screen

## Fitting the statistical model (Optional)

TODO

## Installation

Start by cloning this git repository (likely into your home directory).

ComBind requires access to Glide along with several other Schrodinger tools
and the Schrodinger Python API.

The Schrodinger suite of tools can be accessed on Sherlock by running
`ml chemistry schrodinger`. This will add many of the Schrodinger tools to
your path and sets the SCHRODINGER environmental variable. (Some tools are
not added to your path and you'll need to write out $SCHRODINGER/tool.)
After running this you should be able to run Glide by typing `glide` in the
command line.

You can only access the Schrodinger Python API using their interpretter.
Creating a virtual environment that makes their interpretter the default
python interpretter is the simplest way to do this. To create the environment
and upgrade the relevant packages run the following:

```
cd
$SCHRODINGER/run schrodinger_virtualenv.py schrodinger.ve
source schrodinger.ve/bin/activate
pip install --upgrade numpy sklearn scipy pandas

cd combind
ln -s  ~/schrodinger.ve/bin/activate schrodinger_activate
```

This last line is just there to provide a standardized way to access the
activation script.

Run `source schrodinger_activate` to activate the environment in
the future, you'll need to do this everytime before running ComBind.
This is included in the setup_sherlock script; you can source the
script by running `source setup_sherlock`.


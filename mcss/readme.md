# Maximum Common Substructure Feature

## Overall Methodology

Maximum common substructures (MCSS) are defined as the largest shared set of atoms of the same type present in a pair of ligands. Atom types are defined in the *.typ files in custom_types/.

MCSSs are used to quantify the similarity of two poses by computing the RMSD between
the substructures in each pose. Therefore, ideally, the atom typing and criteria for
when to use the MCSS would be chosen to maximize the possible separation between
native and reference pose pairs.

In addition to being used as a feature, the MCSSs are also used to determine which ligands to used to determine which helper ligands to use. This necessitates that MCSSs be computed for all query-helper ligand pairs. Because of this, the computation is broken into two phases (1) initial MCSS computation using undocked structures and (2) computation of RMSDs between docked poses.

## Code structure

1. The MCSS class represents an individual MCSS and RMSDs. This class handles all of the actual computations and is designed to be a standalone object.

2. The MCSSController class handles decisions about which MCSSs to compute and specifying file paths. This class is therefore specific to the directory structure

## Complications

1. MCSSs are often symetric, which needs to be taken into account in the RMSD computation.

2. Ligand tautomers can change between the undocked and docked poses and amoungst the docked poses.

3. There are lots of pose pairs, making runtimes slow.


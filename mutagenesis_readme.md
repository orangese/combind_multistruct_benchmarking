### Mutagenesis Workflow
#### By: Jason Qin
#### Date: 03/18/19

### Current mutagenesis workflow
1. Download receptor structures for ligands of interest
2. Go through typical Maestro preprocessing (add hydrogens, fix protonation, energy minimize)
3. Decide on reference structure to use for docking of all ligands
4. Generate computational mutants for reference structure
  * Mutagenesis tool in Maestro
  * Rapid torsional scan to choose best rotamer
5. Follow current naming convention to name structures: 
  * [PDB_ID]_[mutation].mae
    * PDB_ID is e.g. 2YDO
    * mutation is e.g. 181S (can have whatever you want as this name, in whatever format you want as long as there aren’t '_' symbols)
6. Align all structures to reference in Maestro
7. Remove ligands from all structures
8. Upload ligands + structures on to Sherlock
9. Prepare ligands + receptors using Combind workflow
10. Dock ligands to receptors
  * Use the “3m” command instead of “3” for prepare_all.py command
11. Run mutant_helpers.py commands to compare docking results across mutants
  * Relevant command examples:
  * $SCHRODINGER/run ~/combind/mutant_helpers.py compare A2AR 2YDO
     * Run this in a batch script
     * Will generate all pairwise RMSDs for all docking results for the 2YDO structure + its mutants 
   * $SCHRODINGER/run ~/combind/mutant_helpers.py aggregate DHFR 3GHW False
     * Will aggregate results for 3GHW structure + its mutants (which mutants are chosen depends on whether control flag is False or True)

### Relevant Directories
* /scratch/PI/rondror/combind/mutagenesis : main directory with mutagenesis experiments
* /scratch/PI/rondror/combind/mutagenesis/[protein]/rmsd : directory containing pairwise pose RMSD results between dockings to different receptors
  * "mutant_helpers.py compare" command will write to this directory 
* /scratch/PI/rondror/combind/mutagenesis/[protein]/rescored : directory containing original pose rankings (i.e. docking ligand to native receptor) and rescored pose rankings rankings (i.e. incorporating docking results from relevant mutants), stored as [ligandName]_rescored.csv
  * "mutant_helpers.py aggregate" command will write to this directory
  * directory also contains the scores of docking for all ligand-receptor pairs, stored in a CSV file for each ligand (e.g. [ligandName]_lig_scores.csv)

### Ideas for future
1. Let the receptor backbone be flexible to better model how the ligand may fit into the binding pocket of the mutant receptors
2. Potential future receptors to characterize
  * http://biosig.unimelb.edu.au/platinum/results
  * Dihydrofolate reductase
     * DH1 - 6 mutant receptors with reasonable Ki
     * MTX - 4 mutant receptors with reasonable Kd
     * http://www.jbc.org/content/270/10/5057.long
     * http://www.jbc.org/content/284/30/20079.long
     * https://pubs.acs.org/doi/10.1021/bi801960h
  * ATP-dependent molecular chaperone HSP82
     * ADP - 5 mutations with reasonable Kd
     * GDM - 4 single mutations with reasonable Kd
     * Dodecin - lots of mutants of same residue
     * LUM
     * RBF
  * Putative ABC transporter amino acid-binding protein
     * 4CS
     * 6CS
  * Thymidylate Synthase
     * DCM
     * UMP
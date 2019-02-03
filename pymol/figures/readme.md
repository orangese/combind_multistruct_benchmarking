
PyMOL scripts used to generate specific figures for the combind manuscript

These are structured as snippets of PyMOL commands to be copied into the window.
See note at beginning of each for instructions.


Figure notes

Figure 1: Using non-structural data to refine ligand binding poses.
 Goal: What problem are we solving?
 Points:
   1. Rescoring: Glide sampling often exceeds scoring.
   2. Ligands often share interactions
   3. We will show how to use this to improve pose scoring.

 Potential display items:
   - Rendering of incorrect glide pose, correct glide pose ranked lower
   - 2D ligand structures
   - Protein binding pocket zoom out
   - Psuedo equation: Physics + Shared Interactions
   - Poses of several crystallographic poses showing shared interactions.
   - Correct Glide pose

 System: AR
   - At least one ligand where ComBind does well.
   - Two well conserved interactions.
   - Stick to hydrogen bonds as these are what we use in the text.
   - Not a GPCR, as we will have two GPCR examples.
   - Preferably a drug target.

Figure 2: Traditional docking algorithms underestimate the rate at which interactions are shared.
 Goal: We show quantitatively that Glide does not adequately capture shared interactions.

 Points:
   1. How do interactions look different in correct v. decoy poses?
   2. How -- roughly -- is similarity calculated?
   3. Show that statistics differ between correct and incorrect poses.

 Potential Display Items:
   - Statistics
   - MCSS overlap
   - Fingerprint

 Decisions:
   - Do we need to show anything here other than the statistics?
      - Pro:
        - Could be unclear what we mean by a maximum common substructure
	- If I was giving a talk I would want to show an example
	- People might not be able to understand without an example
      - Con:
        - The above can be described in the text.
	- It will actually be hard to explain exactly what is going on.

 Previous Feedback:
   - Link between left and right columns are unclear.

 System: DAT, A2AR?

Figure 3: ComBind outperforms standard docking algorithms

 System: None

Figure 4: Case Study: B2AR

Figure 5: Case Study: TRPV1

Figure 6: Prospective application to D2R
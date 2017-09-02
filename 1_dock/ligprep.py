import os
import sys
import numpy as np
from schrodinger.structure import StructureReader
from schrodinger.structutils.analyze import AslLigandSearcher
from schrodinger.structure import StructureWriter
#from schrodinger.infra import mm

def extract_ligands():

    os.system("mkdir -p ligands")
    structs = [f.split('.')[0] for f in os.listdir("raw_maes") if not os.path.exists('ligands/{}_ligand.mae'.format(f.split('.')[0]))]
    #parameters = mm.get_ligand_parameters()	
    #non_ligands = parameters.excluded_residue_names
    #exclude = ['FLF', '056', 'ICO', 'K10', 'YLO', 'YLQ', '4HY', 'MXD', ' T3', 'AV6', ' NK', '2MI', '17W', 'RB1']
    #for ligand in exclude:
	#non_ligands.add(ligand)
    print 'Extracting ligands from {} files'.format(len(structs))
    for s in structs:
        struct = StructureReader("raw_maes/{}.mae".format(s)).next()
        
        asl_searcher = AslLigandSearcher()
        ligands = asl_searcher.search(struct)
        if len(ligands) == 0: 
	   print 'Error: Could not find a ligand for {}'.format(s)
	   out = open("multiple_ligands.out", "a+")
	   out.write("Structure: {}".format(s))
	   out.write(" Number of ligands: 0")
	   out.write("\n")
	   out.close()

	elif len(ligands) == 1:
	   struct.deleteAtoms(ligands[0].atom_indexes)
	   struct.write("raw_maes/{}.mae".format(s))
	   st_writer = StructureWriter("ligands/{}_ligand.mae".format(s))
	   st_writer.append(ligands[0].st)
	   st_writer.close()
	
        else:
 	      for ligand in ligands:
	          if ligands[0].pdbres.strip() != ligand.pdbres.strip():
		     print "Multiple ligands detected in {}. See output file and manually remove them.".format(s)
		     out = open("multiple_ligands.out", "a+")
		     out.write("Structure: {}".format(s))
		     out.write(" Number of ligands: {}".format(len(ligands)))
		     for ligand in ligands:
			out.write(" {}".format(ligand.pdbres))
		     out.write("\n")
		     break
	      else:
		  chain_names = []	  
		  for ligand in ligands:
		      chain_names += [chain._getChainName().strip() for chain in ligand.st.chain]
		  chain_A_ligs = []
		  if "A" in chain_names:
		      for ligand in ligands:
			 if "A" in [chain._getChainName().strip() for chain in ligand.st.chain]:
			    chain_A_ligs.append(ligand)
			 else:
			    struct.deleteAtoms(ligand.atom_indexes)
			    struct.write("raw_maes/{}.mae".format(s)) 	
		      if len(chain_A_ligs) == 1:
		         struct.deleteAtoms(chain_A_ligs[0].atom_indexes)      
                         struct.write("raw_maes/{}.mae".format(s))
                         st_writer = StructureWriter("ligands/{}_ligand.mae".format(s))
                         st_writer.append(chain_A_ligs[0].st)
                         st_writer.close()
	              else:   
		         print "{} has more than 1 ligand in chain A. See output file.".format(s)
		         out = open("multiple_ligands.out", "a+")
	                 out.write("Structure: {}".format(s))
	                 out.write(" Number of ligands: {}".format(len(ligands))) 
	                 for ligand in ligands:
	 	            out.write(" {}".format(ligand.pdbres))
		         out.write("More than 1 ligand in chain A.")
	                 out.write("\n")
	                 out.close()
	          else:	
		      print "{} has no ligand in chain A. See output file.".format(s)
		      out = open("multiple_ligands.out", "a+")
		      out.write("Structure: {}".format(s))
		      out.write(" Number of ligands: {}".format(len(ligands)))
		      for ligand in ligands:
 			 out.write(" {}".format(ligand.pdbres))
	              out.write("No ligand in chain A.")
                      out.write("\n")
		      out.close()

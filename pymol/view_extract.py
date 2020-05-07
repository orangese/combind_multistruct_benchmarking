from pymol import cmd
import sys
import pandas as pd
import glob

iuphar_protein = pd.read_csv('/Users/jpaggi/sherlock/oak/projects/ligand-docking/functional/data/docking/beta2/structures/iuphar_ligands.csv')


iuphar_protein["lig"] = [str(pdb_id)+"_lig" for pdb_id in iuphar_protein["PDB ID"]]
iuphar_protein.set_index("lig")
a = [idx for idx, row in iuphar_protein.iterrows() if re.match(".*gonist.*", row["action"])]
iuphar_protein = iuphar_protein.loc[a]
iuphar_protein["bin_action"] = [0 if "antagonist" in action.lower() else 1 for action in iuphar_protein["action"]]
agonist = iuphar_protein.loc[iuphar_protein["bin_action"]==1]
antagonist = iuphar_protein.loc[iuphar_protein["bin_action"]==0]

m = 10000
structs = glob.glob("*mae")
for s in structs:
    name = s.split('/')[-1]
    name = "".join(name.split('.')[0:-1])
    cmd.load(s, name)
    m -=1
    if m==0:
        break
 
for lig in agonist["lig"]:
    lig_name = lig.split(".")[0]
    cmd.group("agonist combind", "*"+lig_name+"*combind*")
    #cmd.group("agonist top ten", "*"+lig_name+"*ten*")
# cmd.group("agonist combind", "*lig*")
    
for lig in antagonist["lig"]:
    lig_name = lig.split(".")[0]
    cmd.group("antagonist combind", "*"+lig_name+"*combind*")
    #cmd.group("antagonist top ten", "*"+lig_name+"*ten*")

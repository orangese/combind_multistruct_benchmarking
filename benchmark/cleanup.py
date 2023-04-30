import os
import sys
import shutil


def clean(dir_):
    for protein in os.listdir(dir_):
        if not os.path.isdir(protein) or protein == "additional_smiles":
            continue

        ref_ligands = []
        helper_ligands = []

        helpers = os.path.join(dir_, protein, "raw/additional_smiles.smi")
        try:
            with open(helpers, "r") as ligs:
                helper_ligands = [lig.split()[1].rstrip() 
                                  for lig in ligs.readlines()[1:]]
        except FileNotFoundError:
            pass

        for f in os.listdir(os.path.join(dir_, protein)):
            if f.endswith(".out") or f.endswith(".csv"):
                os.remove(os.path.join(dir_, protein, f))
            elif f.startswith("features_"):
                ref_ligands.append(f.strip("features_"))
                shutil.rmtree(os.path.join(dir_, protein, f))

        for todo in ("ligands", "docking"):
            for subdir in os.listdir(os.path.join(dir_, protein, todo)):
                ok = False
                for ref in ref_ligands:
                    if ref in subdir:
                        ok = True
                        break
                for helper in helper_ligands:
                    if helper in subdir:
                        ok = False
                        break
                
                if not ok:
                    shutil.rmtree(os.path.join(dir_, protein, todo, subdir))


if __name__ == "__main__":
    try:
        dir_to_clean = sys.argv[1]
    except IndexError:
        print("args: DIR_TO_CLEAN")
        sys.exit(1)

    clean(dir_to_clean)

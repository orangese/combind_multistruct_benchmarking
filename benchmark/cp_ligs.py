import os
import sys
from shutil import copytree, rmtree

src = sys.argv[1]
dest = sys.argv[2]

for i, prot in enumerate(os.listdir(src)):
    print(f"Doing {i+1} / {len(os.listdir(src))} copying")
    if not os.path.isdir(os.path.join(src, prot)):
        continue
    for lig in os.listdir(os.path.join(src, prot, "ligands")):
        if "lig" not in lig:
            src_p = os.path.join(src, prot, "ligands", lig)
            dest_p = os.path.join(dest, prot, "ligands", lig)
            if os.path.exists(dest_p):
                rmtree(dest_p)
            copytree(src_p, dest_p)

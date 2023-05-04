import os
import shutil
import glob
import sys

if __name__ == "__main__":
    try:
        src = sys.argv[1]
        dest = sys.argv[2]
    except IndexError:
        print("args: SRC DEST")
        sys.exit(1)

    to_copy = ["ligands", "raw", "structures/raw"]

    for p in os.listdir(src):
        try:
            os.makedirs(os.path.join(dest, p))
            os.makedirs(os.path.join(dest, p, "structures"))
        except FileExistsError:
            pass
        fp = os.path.join(src, p)
        for d in to_copy:
            src_d = os.path.join(fp, d)
            # print(src_d, os.path.join(dest, p, d))
            shutil.copytree(src_d, os.path.join(dest, p, d))

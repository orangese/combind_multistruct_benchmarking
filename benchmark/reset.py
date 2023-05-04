import glob
import os
import shutil
import sys

if __name__ == "__main__":
    try:
        src = sys.argv[1]
        dest = sys.argv[2]
        type_ = sys.argv[3]
    except IndexError:
        print("args: SRC DEST TYPE")
        sys.exit(1)

    for f in os.listdir(src):
        try:
            if type_ == "structures/raw":
                src_f = glob.glob(os.path.join(src, f, "structures/raw/AF_*_prot.mae"))[0]
                dest_d = os.path.join(dest, f, "structures/raw/")
                shutil.copy(src_f, dest_d)
            elif type_ == "ligands":
                src_d = os.path.join(src, f, "ligands")
                dest_d = os.path.join(dest, f, "ligands")
                shutil.copytree(src_d, dest_d)
            else:
                raise ValueError("invalid type", type_)
        except IndexError:
            print(f, "failed")


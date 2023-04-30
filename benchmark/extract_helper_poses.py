import shutil
import os
import glob
import sys
import pymol2

def aff_map(map_file):
	mapping = {}
	with open(map_file, "r") as f:
		for t in f.readlines()[1:]:
			try:
				prot, lig = t.split()
			except:
				prot, lig = t.split()[0], None
			mapping[prot] = lig
	return mapping

def highest_aff(p, prot_path, affmap, combind_dir, combind2_dir):
    """gets highest affinity binder by path name. If combind helpers not
    available, gets helpers from second set of ligands (masha's original list)
    """
    prot = prot_path[prot_path.rfind("/") + 1:].rstrip("_prot.mae").strip("AF_")
    # (p is the benchmark, ex NK1R, while prot is the protein, ie 6E59)

    helper = affmap[p]
    if not helper:
        # fetch from backup dir
        combind_dir = combind2_dir
        helper = glob.glob(os.path.join(combind2_dir, p, "docking/*[0-9]-to-AF_*"))[0]
        helper = helper[helper.rfind("/") + 1:]
        helper = helper[:helper.rfind("-to-")]

    dir_ = f"{helper}-to-AF_{prot}"
    lig_path = os.path.join(combind_dir, p, "docking", dir_, dir_ + "_pv.maegz")
    return lig_path, helper


def extract(protlig, out, lig_id, backup=True):
    """protlig should be xxxx-to-AF_xxxx_pv.maegz"""
    if backup:
        shutil.copyfile(out, out.replace("_prot.mae", "_prot.mae.backup"))

    with pymol2.PyMOL() as pymol:
        # load ligand pose
        pymol.cmd.load(protlig, "prot")

        pymol.cmd.copy("lig", f"prot.{lig_id}")
        pymol.cmd.delete("prot")

        # load original structure
        pymol.cmd.load(out)
        pymol.cmd.save(out)
#         input(f"press to coninute, ust did, {protlig} {out} {lig_id}")

if __name__ == "__main__":
    try:
        map_file = sys.argv[1]
        combind = sys.argv[2]
        combind2 = sys.argv[3]
        out = sys.argv[4]
    except IndexError:
        print("args: MAP_FILE COMBIND_DIR1 COMBIND_DIR2 OUT_DIR")
        sys.exit(1)

    affmap = aff_map(map_file)

    for f in os.listdir(out):
        # these should be protein name directories
        prot_to_replace = glob.glob(os.path.join(out, f, "structures/raw/AF_*_prot.mae"))[0]
        # we're going to add the helper ligand pose to this structure

        ref_pose, lig_id = highest_aff(f, prot_to_replace, affmap, combind, combind2)
        extract(ref_pose, prot_to_replace, lig_id)

import glob
import subprocess
import os
import sys
import time

def queue_empty():
    result = subprocess.run("squeue -u $USER", shell=True, capture_output=True)
    out = [s.split() for s in result.stdout.decode().split("\n") if s]
    return len(out) == 1 and out[0][0] == "JOBID"

def wait_until_finish(polltime=60, mock=False):
    if mock:
        return
    time.sleep(10)  # let queue refresh
    while True:
        if queue_empty():
            break
        time.sleep(polltime)
    print("Jobs finished!")

def structprep(pdir, verbose=True, mock=False):
    cmd = "sbatch -p owners --wrap 'combind structprep AF_{}' --output='{}/structprep-%j.out' --job-name='structprep-{}-{}'"
    for p in os.listdir(pdir):
        cwd = os.getcwd()
        fullp = os.path.join(pdir, p)
        os.chdir(fullp)

        for struct in glob.glob(os.path.join("structures/raw/AF_*_prot.mae")):
            struct = struct.rstrip("_prot.mae")
            struct = struct[struct.rfind("_") + 1:]
            tmp_cmd = cmd.format(struct, fullp, struct, p)
            if verbose:
                print(f"Executing {tmp_cmd} in {os.getcwd()}")
            if not mock:
                os.system(tmp_cmd)

        os.chdir(cwd)

    print("Submitted structprep jobs, waiting for completion...")

def ligprep(pdir, verbose=True, mock=False):
    cmd = "sbatch -p owners --wrap 'combind ligprep raw/additional*.smi ligands' --output='{}/ligprep-%j.out' --job-name='ligprep-{}'"
    for p in os.listdir(pdir):
        cwd = os.getcwd()
        fullp = os.path.join(pdir, p)
        os.chdir(fullp)

        tmp_cmd = cmd.format(fullp, p)
        if verbose:
            print(f"Executing {tmp_cmd} in {os.getcwd()}")
        if not mock:
            os.system(tmp_cmd)

        os.chdir(cwd)

    print("Submitted ligprep jobs, waiting for completion...")

def dock(pdir, verbose=True, mock=False):
    cmd = "sbatch -p rondror -t 8:00:00 --wrap 'combind dock docking ligands/*/*.maegz --grid structures/grids/{}/{}.zip' --output='{}/dock-%j.out' --job-name='dock-{}-{}'"
    for p in os.listdir(pdir):
        cwd = os.getcwd()
        fullp = os.path.join(pdir, p)
        os.chdir(fullp)

        for struct in glob.glob(os.path.join("structures/grids/*")):
            struct = struct[struct.rfind("/") + 1:]
            tmp_cmd = cmd.format(struct, struct, fullp, p, struct)
            if verbose:
                print(f"Executing {tmp_cmd} in {os.getcwd()}")
            if not mock:
                os.system(tmp_cmd)

        os.chdir(cwd)

    print("Submitted docking jobs, waiting for completion...")

def featurize(pdir, verbose=True, mock=False):
    cmd = "sbatch -p owners -c 8 --wrap 'combind featurize features_{} docking/*/*{}_pv.maegz --max-poses 25' --output='{}/featurize-%j.out' --job-name='featurize-{}-{}'"

    for p in os.listdir(pdir):
        cwd = os.getcwd()
        fullp = os.path.join(pdir, p)
        os.chdir(fullp)

        for struct in glob.glob(os.path.join("structures/grids/*")):
            struct = struct[struct.rfind("/") + 1:]
            tmp_cmd = cmd.format(struct, struct, fullp, p, struct)
            if verbose:
                print(f"Executing {tmp_cmd} in {os.getcwd()}")
            if not mock:
                os.system(tmp_cmd)

        os.chdir(cwd)

    print("Submitted featurization jobs, waiting for completion...")

def predict(pdir, verbose=True, mock=False):
    cmd = "sbatch -p owners -c 8 --wrap 'combind pose-prediction features_{} out_{}.csv docking/*to-{}/*_pv.maegz' --output='{}/predict-%j.out' --job-name='predict-{}-{}'"

    for p in os.listdir(pdir):
        cwd = os.getcwd()
        fullp = os.path.join(pdir, p)
        os.chdir(fullp)

        for struct in glob.glob(os.path.join("structures/grids/*")):
            struct = struct[struct.rfind("/") + 1:]
            tmp_cmd = cmd.format(struct, struct, struct, fullp, p, struct)
            if verbose:
                print(f"Executing {tmp_cmd} in {os.getcwd()}")
            if not mock:
                os.system(tmp_cmd)

        os.chdir(cwd)

    print("Submitted pose prediction jobs, waiting for completion...")

def rmsd(pdir, verbose=True, mock=False, rmsd=None):
    if rmsd is None:
        rmsd = "/oak/stanford/groups/rondror/projects/ligand-docking/modeled_structures/scripts_debug/compute_rmsd_all.py"
    cmd = "sbatch -p owners --wrap 'python {} {}' --output='{}/rmsd-%j.out' --job-name='rmsd-{}-{}'"
    
    for p in os.listdir(pdir):
        cwd = os.getcwd()
        fullp = os.path.join(pdir, p)
        os.chdir(fullp)

        ref = glob.glob(os.path.join(fullp, "structures", "raw", "AF_*_prot.mae"))[0]
        ref = ref[ref.rfind("/") + 1:]
        ref = ref.rstrip("_prot.mae") 

        structs = []
        for struct in os.listdir(os.path.join(fullp, "structures", "processed")):
            # struct = struct[struct.rfind("/") + 1:]
            # struct = struct.rstrip("_prot.mae")
            if struct != ref:
                structs.append(struct)

        # structs = [struct, struct]  # TODO: fix temp hack
        structs.insert(0, ref)
        tmp_cmd = cmd.format(rmsd, " ".join(structs), fullp, p, ref)
        if verbose:
            print(f"Executing {tmp_cmd} in {os.getcwd()}")
        if not mock:
            os.system(tmp_cmd)

        os.chdir(cwd)

    print("Submitted rmsd jobs, waiting for completion...")


if __name__ == "__main__":
    try:
        pdir = sys.argv[1]
        which = sys.argv[2] if len(sys.argv) > 2 else 'SLDFC'
        mock = bool(sys.argv[-1]) and len(sys.argv) > 3
    except IndexError:
        print("args: BENCHMARK_DIR PIPELINE MOCK_RUN")
        print("PIPELINE: Use S (structprep), L (ligprep), D (dock), F (featurize), P (pose predict)")
        print("to denote which parts of the pipeline to run (no arg to specify all parts)")
        print("MOCK: use third arg to specify that this is a mock run")
        sys.exit(1)

    #if not mock and not queue_empty():
    #    print("slurm queue not empty, exiting")
    #    sys.exit(1)

    print(f"Running combind on benchmark directory {pdir}")
    if mock:
        print("Running in mock mode, no jobs will be submitted")

    if "S" in which:
        print("\n######### combind: structprep ########")
        print("Prepares protein-ligand complex inputs from structures/raw")
        print("Optionally minimizes complexes")
        structprep(pdir, mock=mock)
        wait_until_finish(mock=mock)

    if "L" in which:
        print("\n######### combind: ligprep ###########")
        print("Prepares 3-D structures for ligands to be docked")
        ligprep(pdir, mock=mock)
        wait_until_finish(mock=mock)

    if "D" in which:
        print("\n######### combind: dock ##############")
        print("Docks ligands from prior step to structures from first step")
        dock(pdir, mock=mock)
        wait_until_finish(mock=mock)

    if "F" in which:
        print("\n######### combind: featurize ###########")
        print("Computes features from docked poses for combind statistics")
        featurize(pdir, mock=mock)
        wait_until_finish(mock=mock)

    if "P" in which:
        print("\n######### combind: predict ###########")
        print("Predicts top pose from helper docking interactions")
        predict(pdir, mock=mock)
        wait_until_finish(mock=mock)

    if "R" in which:
        print("\n######### combind: rmsd ###########")
        print("Compute RMSD between predicted poses and crystal ground truth")
        rmsd(pdir, mock=mock)
        wait_until_finish(mock=mock)

    print("Finished!")

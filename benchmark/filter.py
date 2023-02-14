import sys
import os
from timeit import default_timer as timer
import pandas as pd
import pubchempy as pcp
from schrodinger.structure import SmilesStructure

sys.path.insert(1, "../")
from features.mcss import compute_mcss, n_atoms

# units in terms of micromolar for better numerical stability
units = {"µm": 0, "um": 0, "nm": -3, "pm": -6, "fm": -9}


def mcss_size(smiles1, smiles2):
    # Get maximum common substructure size
    m1 = SmilesStructure(smiles1).get3dStructure(require_stereo=False)
    m2 = SmilesStructure(smiles2).get3dStructure(require_stereo=False)

    mcss = compute_mcss(m1, m2, '../features/mcss16.typ')
    mcss_st = SmilesStructure(mcss['st1'][0].split(',')[0].upper()).get2dStructure()
    return n_atoms(mcss_st) / min(n_atoms(m1), n_atoms(m2))


def select_helpers(smiless, n_helpers, chembl_id):
    # smiless is a list of smiles strings sorted by affinity (best first)
    chosen_smiless = set()
    for	smiles in smiless:
        if n_helpers is not None and len(chosen_smiless) == n_helpers:
            return chosen_smiless
        
        try:
            for chosen_smiles in chosen_smiless:
                if mcss_size(smiles, chosen_smiles) > 0.8:
                    break
            else:
                chosen_smiless.add(smiles)
                print(f"({chembl_id}) Total {len(chosen_smiless)} helper ligands so far")
        except ValueError:
            print(f"({chembl_id}) Invalid SMILES '{smiles}', skipping...")

    return chosen_smiless


def filter_from_benchmark(benchmark_dir, n_helper, out_dir):
    def filter_single(f):
        path = os.path.join(benchmark_dir, f)
        data = pd.read_csv(path)
        chembl_id = f.rstrip(".csv")
        outfile = os.path.join(out_dir, "additional_smiles_" + chembl_id + ".smi")

        # Sort by activity
        try:
            data["activity"] = [val * (10 ** units[unit.lower()]) for val, unit 
                                in zip(data["standard_value"], data["standard_units"])]
            data.sort_values("activity", inplace=True)
        except KeyError:
            print(f"No ligand data found for {f}, skipping...")
            with open(outfile, "w+") as dummy:
                dummy.write("SMILES ID\n")
            return

        # Filter top ligands
        ligands = select_helpers(data["canonical_smiles"], n_helper, chembl_id)
        with open(outfile, "w+") as dump:
            dump.write("SMILES ID\n")
            for lig in ligands:
                cid = pcp.get_compounds(lig, "smiles")[0].cid
                dump.write(f"{lig} {cid}\n")

    for f in os.listdir(benchmark_dir):
        filter_single(f)


if __name__ == "__main__":
    try:
        benchmark_dir = sys.argv[1]
        n_helper = int(sys.argv[2])
        out_dir = sys.argv[3]
    except (IndexError, ValueError):
        print("args: BENCHMARK_DIR NUM_HELPER_LIGANDS OUT_DIR")
        sys.exit(1)

    if not os.path.exists(out_dir):
        os.makedirs(out_dir)

    start = timer()
    filter_from_benchmark(benchmark_dir, n_helper, out_dir)
    print(f"Elapsed: {round((timer() - start) / 60, 4)} minutes")

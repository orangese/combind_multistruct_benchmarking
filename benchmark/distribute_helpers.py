import shutil
import sys
import os


def distribute_helpers(smiles_dir, target_dir, benchmark_alias, make_affinity_map):
    with open(benchmark_alias, "r") as benchmark:
        ds = {}
        for line in benchmark.readlines()[1:]:
            name, uniprot, chembl = line.split()
            ds[chembl] = name

    mapping = {}
    for smi_file in os.listdir(smiles_dir):
        # path format: additional_smiles_{CHEMBLID}.smi
        path = os.path.join(smiles_dir, smi_file)
        chembl_id = smi_file.split("_")[2].rstrip(".smi")
        name = ds[chembl_id]

        # target_dir has directories with names of proteins
        if not make_affinity_map:
            shutil.copyfile(path, os.path.join(target_dir, name, "raw/additional_smiles.smi"))
        else:
            with open(path, "r") as f:
                try:
                    info = f.readlines()[1].split()
                    best_id = info[-2]
                except IndexError:
                    best_id = None
                    print(path, "failed")
            mapping[name] = best_id

    if make_affinity_map is not None:
        with open(make_affinity_map, "w+") as d:
            d.write("Prot Lig\n")
            for name, best_id in mapping.items():
                d.write(f"{name} {best_id or ''}\n")

if __name__ == "__main__":
    try:
        smiles_dir = sys.argv[1]
        target_dir = sys.argv[2]
        benchmark_alias = sys.argv[3]
        make_affinity_map = sys.argv[4]
    except IndexError:
        print("args: SMILES_DIR TARGET_DIR BENCHMARK_ALIAS MAKE_AFFINITY_MAP")
        sys.exit(1)

    distribute_helpers(smiles_dir, target_dir, benchmark_alias, make_affinity_map)

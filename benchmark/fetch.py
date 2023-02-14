from timeit import default_timer as timer
import sys
import os

def fetch(benchmark_file, output_dir):
    print(f"Using benchmark from '{benchmark_file}'")
    with open(benchmark_file, "r") as bench:
        benchmark = map(str.rstrip, bench.readlines())

    for i, id_ in enumerate(benchmark):
        print(f"\n({i+1}). SAVING LIGANDS TO {output_dir}/{id_}.csv...")
        os.system(f"python ../chembl/chembl.py query {id_} "
                  f"--output-fname {output_dir}/{id_}.csv")


if __name__ == "__main__":
    try:
        benchmark = sys.argv[1]
        output_dir = sys.argv[2]
    except IndexError:
        print("args: BENCHMARK_FILE OUTPUT_DIR")
        sys.exit(1)

    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    start = timer()
    fetch(benchmark, output_dir)
    print(f"Elapsed: {round((timer() - start) / 60, 4)} minutes")

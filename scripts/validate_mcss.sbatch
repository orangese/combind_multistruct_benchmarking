#!/bin/bash
#SBATCH --partition=rondror
#SBATCH --time=02:00:00
#SBATCH --tasks=4 --cpus-per-task=1 --ntasks-per-socket=2
#SBATCH --array=0-9%3
#SBATCH --open-mode=append
#SBATCH --mem=8GB
#SBATCH --out=/home/users/jpaggi/validate.out
#SBATCH --err=/home/users/jpaggi/validate.err

mode=$1

. setup_combind

data=/scratch/PI/rondror/combind/bpp_data
end=$((4*SLURM_ARRAY_TASK_ID+4))
for prot in $(ls $data | head -n $end | tail -n 4); do
    $SCHRODINGER/run $COMBINDHOME/1_dock/prepare_all.py 'v' $prot &
done
wait

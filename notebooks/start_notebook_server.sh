#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=jupyter
#SBATCH --tasks=1 --cpus-per-task=1
#SBATCH --ntasks-per-socket=1
##SBATCH -N 1 -n 8
##SBATCH --qos=rondror_high --partition=rondror
#SBATCH --partition=rondror
#SBATCH --output=server_host.out
#SBATCH --dependency=singleton

module load schrodinger/2017-3
export JUPYTER_RUNTIME_DIR=$SCRATCH

hostname
jupyter-notebook --port=8890 --no-browser


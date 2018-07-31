#!/bin/bash
#SBATCH --time=08:00:00
#SBATCH --mem=8GB
#SBATCH --job-name=jupyter
#SBATCH --tasks=1 --cpus-per-task=1
#SBATCH --ntasks-per-socket=1
#SBATCH --partition=rondror
#SBATCH --output=server_host.out
#SBATCH --dependency=singleton

jupyter-notebook --port=8890 --no-browser


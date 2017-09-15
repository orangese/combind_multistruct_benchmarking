#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --mem=16GB
#SBATCH --job-name=jupyter
#SBATCH --tasks=1
#SBATCH --ntasks-per-socket=1
##SBATCH -N 1 -n 8
##SBATCH --qos=rondror_high --partition=rondror
#SBATCH --partition=rondror
#SBATCH --output=server_host.out
#SBATCH --dependency=singleton

hostname
jupyter-notebook --port=8895 --no-browser


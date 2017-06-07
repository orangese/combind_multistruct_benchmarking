#!/bin/bash
#SBATCH --time=8:00:00
#SBATCH --mem=128GB
#SBATCH --job-name=jupyter
#SBATCH --tasks=4
##SBATCH -N 1 -n 8
##SBATCH --qos=rondror_high --partition=rondror
#SBATCH --partition=rondror
#SBATCH --output=server_host.out
#SBATCH --dependency=singleton

hostname
jupyter-notebook --port=8893 --no-browser


#!/bin/bash

#SBATCH --job-name=V5
#SBATCH --qos=debug
#SBATCH --time=00:05:00
#SBATCH --constraint=cpu
#SBATCH --nodes=1
#SBATCH --account=m3972

cd $SLURM_SUBMIT_DIR
module load vasp/6.4.3-cpu

./twoLevelSystem.sh


#!/bin/bash
#SBATCH -Jlymph_final
#SBATCH -Agts-aalexeev3-afrench
#SBATCH -N1
#SBATCH --cpus-per-task=16
#SBATCH -t3-0	# formats: "minutes""minutes:seconds""hours:minutes:seconds""days-hours""days-hours:minutes""days-hours:minutes:seconds".
#SBATCH -qinferno	# embers and inferno
#SBATCH -orun.stdout
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=apoorghani3@gatech.edu

export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK
./lbls
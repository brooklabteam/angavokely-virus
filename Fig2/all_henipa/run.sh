#!/bin/bash
#SBATCH --job-name=henipatree
#SBATCH --account=fc_sirmodel
#SBATCH --partition=savio3
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=7
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

module load vim/7.4
module load emacs/25.1
module load python/3.6
module load java/1.8.0_121
module load cmake/3.15.1
module load gcc/7.4.0
module load openmpi/4.0.1-gcc

/global/home/users/cbrook/raxml-ng/raxml-ng-mpi --all --msa all_henipa_align.fasta --model GTR+I+G4 --prefix T3  --seed 1 --threads 7 --bs-metric fbp,tbe

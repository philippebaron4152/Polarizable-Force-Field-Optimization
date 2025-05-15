#!/bin/bash
#SBATCH --job-name=solution
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=02:00:00
#SBATCH --mail-user=pb4152@princeton.edu
#SBATCH --mail-type=end
#SBATCH --mail-type=fail

module purge
module load intel-oneapi/2024.2
module load intel-mpi/oneapi/2021.13
module load intel-mkl/2024.2
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK


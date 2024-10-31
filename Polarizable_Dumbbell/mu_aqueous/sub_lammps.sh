#!/bin/bash
#SBATCH --job-name=solution
#SBATCH --nodes=1
#SBATCH --ntasks=16
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --time=12:00:00
#SBATCH --mail-user=pb4152@princeton.edu
#SBATCH --mail-type=end
#SBATCH --mail-type=fail
#SBATCH --constraint=haswell

module purge
module load intel/19.1.1.217
module load intel-mpi/intel/2019.7
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

srun $HOME/.local/bin/lmp_della -in in.lammps

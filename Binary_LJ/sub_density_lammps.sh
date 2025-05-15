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

module load anaconda3/2023.3
conda activate myenv

path=/scratch/gpfs/pb4152/Example/Polarizable-Force-Field-Optimization/Binary_LJ/

srun $HOME/.local/bin/lmp_della -in in_density.lammps
wait
python $path/compute_density.py
rho=$(head -n 1 density.txt)
N=N_ATOMS

cd ../fep
python $path/crystal_setup.py -p "$PWD" -n ${N} -r ${rho} -l "FCC" -b ""
sbatch sub_solid_lammps.sh
cd -

cd ../lattice
python $path/crystal_setup.py -p "$PWD" -n ${N} -r ${rho} -l "FCC" -b ""
sbatch sub_solid_lammps.sh
cd -

n=1
while read line; do 
    p=../springs/k${n}/
    cd ${p}
    python $path/crystal_setup.py -p "$PWD" -n ${N} -r ${rho} -l "FCC" -b ""
    sbatch sub_solid_lammps.sh
    cd -
    n=$(($n+1))
done < $path/lambdas_solid.txt
#!/bin/bash

base_path=$PWD
n_lambdas=$(wc -l < "lambdas.txt")

AA_params=($(sed -n '1p' parameters.txt))
AB_params=($(sed -n '2p' parameters.txt))
BB_params=($(sed -n '3p' parameters.txt))

eps_AA=${AA_params[0]}
sigma_AA=${AA_params[1]}
eps_AB=${AB_params[0]}
sigma_AB=${AB_params[1]}
eps_BB=${BB_params[0]}
sigma_BB=${BB_params[1]}

START=10
STOP=90
INCREMENT=20

for (( N_A=${START}; N_A<=${STOP}; N_A=$(($N_A+$INCREMENT)) )); do
    n=1
    while read line; do
        path=$base_path/aqueous/${N_A}A/lambda_${n}/
        if [ ! -d $path ]; then 
            mkdir -p $path
        fi

        cp ./{in_fluid.lammps,sub_fluid_lammps.sh,*.xyz,mixture.inp} $path
        cd $path
        sed -i "s/N_A/${N_A}/g" mixture.inp
        $base_path/packmol < mixture.inp > /dev/null
        python $base_path/packmol_to_lammps.py -i binary_lj.xyz -o data.lammps -m mixture.inp -b "" -a "" 
        sed -i "s/LAMBDA/${line}/g" in_fluid.lammps
        sed -i "s/EPS_AA/${eps_AA}/g" in_fluid.lammps
        sed -i "s/EPS_AB/${eps_AB}/g" in_fluid.lammps
        sed -i "s/EPS_BB/${eps_BB}/g" in_fluid.lammps
        sed -i "s/SIGMA_AA/${sigma_AA}/g" in_fluid.lammps
        sed -i "s/SIGMA_AB/${sigma_AB}/g" in_fluid.lammps
        sed -i "s/SIGMA_BB/${sigma_BB}/g" in_fluid.lammps
        cd - > /dev/null

        touch $path/in_.lammps
        cp ./in_fluid_rerun.lammps $path/in_.lammps
        sed -i "s/LAMBDA/${line}/g" $path/in_.lammps
        sed -i "s/LOG/log_.lammps/g" $path/in_.lammps

        sed -i "s/EPS_AA/${eps_AA}/g" $path/in_.lammps
        sed -i "s/EPS_AB/${eps_AB}/g" $path/in_.lammps
        sed -i "s/EPS_BB/${eps_BB}/g" $path/in_.lammps
        sed -i "s/SIGMA_AA/${sigma_AA}/g" $path/in_.lammps
        sed -i "s/SIGMA_AB/${sigma_AB}/g" $path/in_.lammps
        sed -i "s/SIGMA_BB/${sigma_BB}/g" $path/in_.lammps

        sed -i '$a\wait' $path/sub_fluid_lammps.sh
        sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_.lammps' $path/sub_fluid_lammps.sh

        if [ $n -lt $n_lambdas ]; then
            touch $path/in_up.lammps
            line_number=$(($n+1))
            lambda_up=$(sed -n "${line_number}p" "lambdas.txt")
            cp ./in_fluid_rerun.lammps $path/in_up.lammps
            sed -i "s/LAMBDA/${lambda_up}/g" $path/in_up.lammps
            sed -i "s/LOG/log_up.lammps/g" $path/in_up.lammps 

            sed -i "s/EPS_AA/${eps_AA}/g" $path/in_up.lammps
            sed -i "s/EPS_AB/${eps_AB}/g" $path/in_up.lammps
            sed -i "s/EPS_BB/${eps_BB}/g" $path/in_up.lammps
            sed -i "s/SIGMA_AA/${sigma_AA}/g" $path/in_up.lammps
            sed -i "s/SIGMA_AB/${sigma_AB}/g" $path/in_up.lammps
            sed -i "s/SIGMA_BB/${sigma_BB}/g" $path/in_up.lammps

            sed -i '$a\wait' $path/sub_fluid_lammps.sh
            sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_up.lammps' $path/sub_fluid_lammps.sh
        fi

        if [ $n -gt 1 ]; then
            touch $path/in_down.lammps
            line_number=$(($n-1))
            lambda_down=$(sed -n "${line_number}p" "lambdas.txt")
            cp ./in_fluid_rerun.lammps $path/in_down.lammps
            sed -i "s/LAMBDA/${lambda_down}/g" $path/in_down.lammps
            sed -i "s/LOG/log_down.lammps/g" $path/in_down.lammps

            sed -i "s/EPS_AA/${eps_AA}/g" $path/in_down.lammps
            sed -i "s/EPS_AB/${eps_AB}/g" $path/in_down.lammps
            sed -i "s/EPS_BB/${eps_BB}/g" $path/in_down.lammps
            sed -i "s/SIGMA_AA/${sigma_AA}/g" $path/in_down.lammps
            sed -i "s/SIGMA_AB/${sigma_AB}/g" $path/in_down.lammps
            sed -i "s/SIGMA_BB/${sigma_BB}/g" $path/in_down.lammps

            sed -i '$a\wait' $path/sub_fluid_lammps.sh
            sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_down.lammps' $path/sub_fluid_lammps.sh 
        fi
        
        echo "LAMBDA = " ${line} ", SETUP DONE"

        cd $path
        sbatch sub_fluid_lammps.sh
        cd - > /dev/null

        n=$(($n+1))
    done < "lambdas.txt"
done

for N in 4 6; do
    n_atoms=$(((2*${N})*(2*${N})*(2*${N})/2))
    path_density=$base_path/crystal/${n_atoms}/density/

    if [ ! -d $path_density ]; then 
        mkdir -p $path_density
    fi

    cp ./{in_density.lammps,sub_density_lammps.sh} ${path_density}
    sed -i "s/N_ATOMS/${N}/g" $path_density/sub_density_lammps.sh
    sed -i "s/EPS_AA/${eps_AA}/g" $path_density/in_density.lammps
    sed -i "s/SIGMA_AA/${sigma_AA}/g" $path_density/in_density.lammps
    cd $path_density
    python $base_path/crystal_setup.py -p "$PWD" -n ${N} -r 1500 -l "FCC" -b "" 
    sbatch sub_density_lammps.sh
    cd -

    path_fep=$base_path/crystal/${n_atoms}/fep/
    path_lattice=$base_path/crystal/${n_atoms}/lattice/

    if [ ! -d $path_fep ]; then 
        mkdir -p $path_fep
    fi

    cp ./{in_solid.lammps,in_solid_rerun.lammps,sub_solid_lammps.sh} ${path_fep}
    sed -i "s/LAMBDA/0.0/g" $path_fep/in_solid.lammps
    sed -i "s/SPRING_SCALING/1.0/g" $path_fep/in_solid.lammps
    sed -i "s/LAMBDA/1.0/g" $path_fep/in_solid_rerun.lammps

    sed -i "s/EPS_AA/${eps_AA}/g" $path_fep/in_solid.lammps
    sed -i "s/SIGMA_AA/${sigma_AA}/g" $path_fep/in_solid.lammps
    sed -i "s/EPS_AA/${eps_AA}/g" $path_fep/in_solid_rerun.lammps
    sed -i "s/SIGMA_AA/${sigma_AA}/g" $path_fep/in_solid_rerun.lammps

    sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_solid.lammps' $path_fep/sub_solid_lammps.sh 
    sed -i '$a\wait' $path_fep/sub_solid_lammps.sh 
    sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_solid_rerun.lammps' $path_fep/sub_solid_lammps.sh 

    if [ ! -d $path_lattice ]; then 
        mkdir -p $path_lattice
    fi
    cp ./{in_lattice.lammps,sub_solid_lammps.sh} ${path_lattice}
    sed -i "s/EPS_AA/${eps_AA}/g" $path_lattice/in_lattice.lammps
    sed -i "s/SIGMA_AA/${sigma_AA}/g" $path_lattice/in_lattice.lammps
    sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_lattice.lammps' $path_lattice/sub_solid_lammps.sh 

    n=1
    while read line; do 
        path_springs=$base_path/crystal/${n_atoms}/springs/k${n}/
        if [ ! -d $path_springs ]; then 
            mkdir -p $path_springs
        fi
        cp ./{in_solid.lammps,sub_solid_lammps.sh} ${path_springs}
        sed -i "s/LAMBDA/1.0/g" $path_springs/in_solid.lammps
        sed -i "s/EPS_AA/${eps_AA}/g" $path_springs/in_solid.lammps
        sed -i "s/SIGMA_AA/${sigma_AA}/g" $path_springs/in_solid.lammps
        sed -i "s/SPRING_SCALING/${line}/g" $path_springs/in_solid.lammps
        sed -i '$a\srun $HOME/.local/bin/lmp_della -in in_solid.lammps' $path_springs/sub_solid_lammps.sh 
        n=$(($n+1))
    done < "lambdas_solid.txt"
done
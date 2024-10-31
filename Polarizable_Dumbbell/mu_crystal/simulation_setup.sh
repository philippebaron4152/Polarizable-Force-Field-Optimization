#!/bin/bash

N_pair=6
N=$(((2*${N_pair})*(2*${N_pair})*(2*${N_pair})/2))

run_num=1
N_rep=1
rerun=false

fep=true

for (( i=1; i<=${N_rep}; i++ ))
do
	D=${N}/${i}/

	if [ ! -d "${D}" ]; then
                mkdir -p ${D}
        fi
	
	cd ${D}
	python /scratch/gpfs/pb4152/Scripts/crystal_setup.py -p "$PWD" -n ${N_pair} -r 2490 -l "CP1" -b "A A1" -d "A A1" -sy "" -m "A" --polarizable
	cd -

	if ${fep}; then
		if ! ${rerun}; then
			DIR=${D}/fep/run${run_num}
			DIR2=${D}/fep/lattice	

			if [ ! -d "${DIR}" ]; then
			      mkdir -p ${DIR}
			fi
			
			if [ ! -d "${DIR2}" ]; then
                              mkdir -p ${DIR2}
                        fi

			SEED=$((27856+${i}*${i}))

			cp ./{in_solid.lammps,sub_solid_lammps.sh} ${DIR}
			cp ${D}/data_crystal.lammps ${DIR}

			sed -i "s/LAMBDA/0.0/g" ${DIR}/in_solid.lammps
			sed -i "s/PHI/0.0/g" ${DIR}/in_solid.lammps
			sed -i "s/SPRING_SCALING/1.0/g" ${DIR}/in_solid.lammps
			sed -i "s/REPLACE_SEED/${SEED}/g" ${DIR}/in_solid.lammps
			sed -i "s/KSPACE_STYLE//g" ${DIR}/in_solid.lammps
			sed -i "s#PAIR_STYLE#lj/cut#g" ${DIR}/in_solid.lammps

			cp ./{in_lattice.lammps,sub_lattice_lammps.sh} ${DIR2}
			cp ${D}/data_crystal.lammps ${DIR2}
			cd $DIR2
			sbatch sub_lattice_lammps.sh
			cd -
		else
			DIR=${D}/fep/rerun
                 
			if [ ! -d "${DIR}" ]; then
					mkdir -p ${DIR}
			fi
		
			SEED=$((27856+${i}*${i}))

			cp ./{in_rerun.lammps,sub_rerun_lammps.sh} ${DIR}
			cp ${D}/data_crystal.lammps ${DIR}
			cp ${D}/fep/run${run_num}/dump.coord ${DIR}

			sed -i "s/LAMBDA/1.0/g" ${DIR}/in_rerun.lammps
			sed -i "s/PHI/1.0/g" ${DIR}/in_rerun.lammps
			sed -i "s/SPRING_SCALING/1.0/g" ${DIR}/in_rerun.lammps
			sed -i "s/REPLACE_SEED/${SEED}/g" ${DIR}/in_rerun.lammps
		fi
	else
		n=1
		while read line; do
			DIR=${D}/springs/x${n}/

			if [ ! -d "${DIR}" ]; then
					mkdir -p ${DIR}
			fi

			cp ./{in_solid.lammps,sub_solid_lammps.sh} ${DIR}
			cp ${D}/data_crystal.lammps ${DIR}

			SEED=$((34567+${n}*${n}+${i}))

			sed -i "s/LAMBDA/1.0/g" ${DIR}/in_solid.lammps
			sed -i "s/PHI/1.0/g" ${DIR}/in_solid.lammps
			sed -i "s/SPRING_SCALING/${line}/g" ${DIR}/in_solid.lammps
			sed -i "s/REPLACE_SEED/${SEED}/g" ${DIR}/in_solid.lammps
			sed -i "s/KSPACE_STYLE/kspace_style 	pppm 1.0e-6/g" ${DIR}/in_solid.lammps
			sed -i "s#PAIR_STYLE#lj/cut/coul/long#g" ${DIR}/in_solid.lammps

			n=$(($n+1))             
		done < "lambdas.txt"
	fi
done

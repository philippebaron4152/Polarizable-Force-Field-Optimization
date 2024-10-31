#!/bin/bash
N_B=1350

vdw=true
run_num=1
prev_run=$((${run_num}-1))

SCRIPTS=/scratch/gpfs/pb4152/Scripts

if [ ${run_num} -eq 1 ]; then
        IN_FILE=data.lammps
        OUT_FILE=out.lammps
elif [ ${run_num} -eq 2 ]; then
        IN_FILE=out.lammps
        OUT_FILE=run2_out.lammps
else
        IN_FILE=run${prev_run}_out.lammps
        OUT_FILE=run${run_num}_out.lammps
fi

START=30
STOP=30

for (( N_A=${START}; N_A<=${STOP}; N_A=$((${N_A}+20)) )) 
do
	D=${N_B}B/${N_A}A/
	
	if [ ! -d "${D}" ]; then
		mkdir -p ${D}
	fi

	if ${vdw}; then
		n=1
		while read line; do
			DIR=${D}vdw/k_${n}/run${run_num}
			if [ ! -d "${DIR}" ]; then
				mkdir -p ${DIR}
			fi

			if [ ${run_num} -eq 1 ]; then 
				cp ./{mixture.inp,*.xyz} ${DIR}
				sed -i "s/N_A/${N_A}/g" ${DIR}/mixture.inp
				sed -i "s/N_B/${N_B}/g" ${DIR}/mixture.inp
				cd ${DIR}
				 ~/software/packmol-20.14.0/packmol < mixture.inp
				python $SCRIPTS/packmol_to_lammps_polarizable.py -i binary_lj.xyz -o data.lammps -m mixture.inp  -a "" -d "A A1 C C1" -s "" -b "A A1 C C1" -mol "A B C"
				cd -
			else
				cp ${D}vdw/k_${n}/run${prev_run}/${IN_FILE} ${DIR}
			fi
				
			echo "$PWD"	
			cp ./{sub_lammps.sh,in.lammps} ${DIR}
			sed -i "s/LAMBDA/${line}/g" ${DIR}/in.lammps
			sed -i "s/PHI/0.0/g" ${DIR}/in.lammps
			sed -i "s/INFILE/${IN_FILE}/g" ${DIR}/in.lammps
			sed -i "s/OUTFILE/${OUT_FILE}/g" ${DIR}/in.lammps

			n=$((${n}+1))
		done < "lambdas_vdw.txt"
	else
		n=1
		while read line; do
			DIR=${D}coul/k_${n}/run${run_num}
			if [ ! -d "${DIR}" ]; then
				mkdir -p ${DIR}
			fi

			if [ ${run_num} -eq 1 ]; then
				cp ./{mixture.inp,*.xyz} ${DIR}
				sed -i "s/N_A/${N_A}/g" ${DIR}/mixture.inp
				sed -i "s/N_B/${N_B}/g" ${DIR}/mixture.inp
				cd ${DIR}
					~/software/packmol-20.14.0/packmol < mixture.inp
				python $SCRIPTS/packmol_to_lammps_polarizable.py -i binary_lj.xyz -o data.lammps -m mixture.inp  -a "" -d "A A1 C C1" -s "" -b "A A1 C C1" -mol "A B C"
				cd -
			else
				cp DIR=${D}coul/k_${n}/run${prev_run}/${INFILE} ${DIR}
			fi

			echo "$PWD"     
			cp ./{sub_lammps.sh,in.lammps} ${DIR}
			sed -i "s/LAMBDA/1.0/g" ${DIR}/in.lammps
			sed -i "s/PHI/${line}/g" ${DIR}/in.lammps
			sed -i "s/INFILE/${IN_FILE}/g" ${DIR}/in.lammps
			sed -i "s/OUTFILE/${OUT_FILE}/g" ${DIR}/in.lammps

			n=$((${n}+1))
		done < "lambdas_coul.txt"
	fi

done


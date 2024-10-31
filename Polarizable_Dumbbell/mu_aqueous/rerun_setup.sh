#!/bin/bash
N_B=1350

vdw=false
run_num=1
prev_run=$((${run_num}-1))

SCRIPTS=/scratch/gpfs/pb4152/Scripts

START=30
STOP=30
INCREMENT=0

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

for (( N_A=${START}; N_A<=${STOP}; N_A=$((${N_A}+${INCREMENT})) )) 
do
	D=${N_B}B/${N_A}A/
	
	if ${vdw}; then
		n=1
		while read line; do
			DIR=${D}vdw/k_${n}/reruns

			if [ ! -d "${DIR}" ]; then
				mkdir -p ${DIR}
			fi

			cp ${D}vdw/k_${n}/run${run_num}/{${IN_FILE},dump.coord} ${DIR}
			
			n2=1
			while read line_; do

				if [ ${n2} -ge $((${n}-1)) && ${n2} -le $((${n}+1)) ]; then
					infile=${DIR}/in_k_${n2}.lammps
					touch ${infile}
					cp in_rerun.lammps ${infile}
					sed -i "s/LAMBDA/${line_}/g" ${infile}
					sed -i "s/PHI/0.0/g" ${infile}
					sed -i "s/LOG/k_${n2}_log.lammps/g" ${infile}
					sed -i "s/INFILE/${IN_FILE}/g" ${infile}

					touch ${DIR}/sub_k_${n2}.sh
					cp sub_rerun_lammps.sh ${DIR}/sub_k_${n2}.sh
					sed -i "s/REPLACE/in_k_${n2}.lammps/g" ${DIR}/sub_k_${n2}.sh
				fi

				n2=$((${n2}+1))
			done < "lambdas_vdw.txt"
			n=$((${n}+1))
		done < "lambdas_vdw.txt"
	else
		n=1
		while read line; do
			DIR=${D}coul/k_${n}/reruns

			if [ ! -d "${DIR}" ]; then
				mkdir -p ${DIR}
			fi

			cp ${D}coul/k_${n}/run${run_num}/{${IN_FILE},dump.coord} ${DIR}
			
			n2=1
			while read line_; do

				if [ ${n2} -ge $((${n}-1)) && ${n2} -le $((${n}+1)) ]; then
					infile=${DIR}/in_k_${n2}.lammps
					touch ${infile}
					cp in_rerun.lammps ${infile}
					sed -i "s/LAMBDA/1.0/g" ${infile}
					sed -i "s/PHI/${line_}/g" ${infile}
					sed -i "s/LOG/k_${n2}_log.lammps/g" ${infile}
					sed -i "s/INFILE/${IN_FILE}/g" ${infile}

					touch ${DIR}/sub_k_${n2}.sh
					cp sub_rerun_lammps.sh ${DIR}/sub_k_${n2}.sh
					sed -i "s/REPLACE/in_k_${n2}.lammps/g" ${DIR}/sub_k_${n2}.sh
				fi

				n2=$((${n2}+1))
			done < "lambdas_coul.txt"
			n=$((${n}+1))
		done < "lambdas_coul.txt"
	fi

done


#!/bin/bash

N_pair=6
N=$(((2*${N_pair})*(2*${N_pair})*(2*${N_pair})/2))

N_rep=1
run_num=1

fep=false
rerun=false
ti=true

for (( i=1; i<=${N_rep}; i++ ))
do
	D=${N}/${i}/

	if ${fep}; then
		if $rerun; then 
			DIR=${D}/fep/rerun
			cd $DIR
               		sbatch sub_rerun_lammps.sh
                	cd -
		else
			DIR=${D}/fep/run${run_num}
                        cd $DIR
                        sbatch sub_solid_lammps.sh
                        cd -
		fi
	elif ${ti}; then 
		n=1
                while read line; do
                        DIR=${D}/ti/x${n}/
                        cd $DIR
                        sbatch sub_solid_lammps.sh
                        cd -
                        n=$(($n+1))
                done < "lambdas_coul.txt"
	else
		n=1
                while read line; do
                        DIR=${D}/springs/x${n}/
			cd $DIR
			sbatch sub_solid_lammps.sh
			cd -
                        n=$(($n+1))             
                done < "lambdas.txt"
	fi
done

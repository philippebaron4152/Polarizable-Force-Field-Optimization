#!/bin/bash
N_B=1350
N_A=30

vdw=false
run_num=1
rerun=true

D=${N_B}B/${N_A}A/

if $rerun; then
	if ${vdw}; then
                n=1
                while read line; do
                        DIR=${D}vdw/k_${n}/reruns
                        cd $DIR
			n2=1
			while read line_; do 
				sbatch sub_k_${n2}.sh
				n2=$((${n2}+1))
			done < ../../../../../"lambdas_vdw.txt"
                        cd -
                        n=$((${n}+1))
                done < "lambdas_vdw.txt"
        else
                n=1
                while read line; do
                        DIR=${D}coul/k_${n}/reruns
                        cd $DIR
                        sbatch sub_k_${n}.sh
                        cd -
                        n=$((${n}+1))
                done < "lambdas_coul.txt"
        fi
else
	if ${vdw}; then
		n=1
		while read line; do
			DIR=${D}vdw/k_${n}/run${run_num}
			cd $DIR
			sbatch sub_lammps.sh
			cd -
			n=$((${n}+1))
		done < "lambdas_vdw.txt"
	else
		n=1
		while read line; do
			DIR=${D}coul/k_${n}/run${run_num}
			cd $DIR
			sbatch sub_lammps.sh
			cd -
			n=$((${n}+1))
		done < "lambdas_coul.txt"
	fi
fi


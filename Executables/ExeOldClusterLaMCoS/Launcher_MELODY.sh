#!/bin/sh
#PBS -N MELODY_MySim
#PBS -l nodes=1:ppn=20,mem=4GB,walltime=96:00:00
#PBS -j oe
#PBS -M surname.name@insa-lyon.fr
#PBS -m abe
#PBS -o LOG_MELODY.log
#PBS -V

module load intel/11.1.080
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=$PBS_NUM_PPN
./MELODY_2D 0

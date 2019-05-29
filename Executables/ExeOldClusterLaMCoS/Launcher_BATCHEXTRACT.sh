#!/bin/sh
#PBS -N BATCHEXTRACT_MySim
#PBS -l nodes=1:ppn=1,mem=4GB,walltime=4:00:00
#PBS -j oe
#PBS -M surname.name@insa-lyon.fr
#PBS -m abe
#PBS -o LOG_BATCHEXTRACT.log
#PBS -V

module load intel/11.1.080
cd $PBS_O_WORKDIR
export OMP_NUM_THREADS=$PBS_NUM_PPN
./BATCHEXTRACT 0 1 200

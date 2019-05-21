#!/bin/bash
#SBATCH -J So016
#SBATCH --partition=parallel
#SBATCH --time='48:00:00'
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=16
#SBATCH --mem-per-cpu=4000
#SBATCH --mail-type=end
#SBATCH --mail-user=guilhem.mollon@insa-lyon.fr
#SBATCH -o output_%j.txt
#SBATCH -e stderr_%j.txt

echo "Starting at `date`"
echo "Running on hosts : $SLURM_NODELIST"
echo "Running on $SLURM_NNODES nodes."
echo "Running on $SLURM_NPROCS processors."

export OMP_NUM_THREAD=16
module load gcc/4.9.2
./MELODY_2D 0
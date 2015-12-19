#!/bin/bash
#SBATCH --error=OpenMP_SM_Regress.err
#SBATCH --output=OpenMP_SM_Regress.out
#SBATCH --job-name=CSE603-OpenMP_SM_Regress
#SBATCH --mail-user=u2@buffalo.edu
#SBATCH --mail-type=ALL
#SBATCH --time=1:00:00 # 24 hours to run this 
##SBATCH --cpus-per-task=1
#SBATCH --nodes=2
#SBATCH --partition=debug
#SBATCH --mem=48000
#SBATCH --tasks-per-node=8
##SBATCH --constraint=CPU-E5-2660















tic=`date +%s`
echo "Start Time = "`date`

#construct nodefile
SLURM_NODEFILE=my_slurm_nodes.$$
srun hostname | sort > $SLURM_NODEFILE

# load modules
module load intel/16.0
module load cilk/5.4.6
# cd to directory from which job was submitted
cd $SLURM_SUBMIT_DIR

# compute number of processors
NPROCS=`cat $SLURM_NODEFILE | wc -l`

# set intel-mpi environment variables
# turn debugging up a bit
P=16 ## number of processors
_Q=1000 
_N=300
_K=6
runMax=20;
export OMP_NUM_THREADS=$P


# launch mainDriver.cpp
# NOTE: This requires you be in the workdir
echo "Launching .."
module list 
##srun -n $NPROCS mainDriver.cpp
srun  icc -openmp -o3 -std=c++11 mainDriver.cpp solution.cpp dataPreparation.cpp -o test1
srun ./test1 1 1 1 20
echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds"

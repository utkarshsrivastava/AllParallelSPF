#!/bin/bash
#SBATCH --error=OpenMP_SM_Regress.err
#SBATCH --output=OpenMP_SM_Regress.out
#SBATCH --job-name=CSE603-OpenMP_SM_Regress
#SBATCH --mail-user=afarasat@buffalo.edu
#SBATCH --mail-type=END
#SBATCH --time=24:00:00 # 24 hours to run this 
##SBATCH --cpus-per-task=2
#SBATCH --nodes=1 
##SBATCH --partition=debug
#SBATCH --mem=8000
#SBATCH --tasks-per-node=12















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
P=1 ## number of processors
_Q=200000 
_N=5000
_K=8
runMax=20;
export OMP_NUM_THREADS=$P


# launch mainDriver.cpp
# NOTE: This requires you be in the workdir
echo "Launching .."
module list 
##srun -n $NPROCS mainDriver.cpp
srun  icc -openmp -o3 -std=c++11 mainDriver.cpp solution.cpp dataPreparation.cpp -o test6
srun ./test6 $_Q $_N $_K $P $runMax
echo "All Done!"

echo "End Time = "`date`
toc=`date +%s`

elapsedTime=`expr $toc - $tic`
echo "Elapsed Time = $elapsedTime seconds"

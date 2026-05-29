#!/bin/bash
#SBATCH --job-name=Geq_D3D
#SBATCH --account OPEN-33-33
#SBATCH --partition=qcpu
#SBATCH --nodes=16
#SBATCH --ntasks=576
#SBATCH --cpus-per-task=1
#SBATCH --time=48:00:00

echo "Starting job $SLURM_JOB_ID on $SLURM_JOB_NODELIST"
module load MCR/R2021a

##  change to scratch directory
# SCR=/scratch/temp/$USER
SLURM_SUBMIT_DIR=/ramdisk/$SLURM_JOB_ID 
mkdir -p $SLURM_SUBMIT_DIR ; 
echo ramdisk folder created!

# run all 144 tasks in parallel with a single srun
srun --ntasks=576 --cpus-per-task=1 ./wrapper_203301_30L_1DF06.sh

echo "All tasks finished."

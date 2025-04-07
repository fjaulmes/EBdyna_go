#!/bin/bash
#SBATCH --job-name Geq_D3D
#SBATCH --account OPEN-33-60
#SBATCH --partition qcpu
#SBATCH --nodes 4
#SBATCH --ntasks-per-node 36
#SBATCH --ntasks 144
#SBATCH --cpus-per-task=1
#SBATCH --time 16:00:00
                         # 16 hours wall-clock 
                         # time allowed for this job
                         # 8 nodes for this job
module load MCR/R2021a


##  change to scratch directory
# SCR=/scratch/temp/$USER
SLURM_SUBMIT_DIR=/ramdisk/$SLURM_JOB_ID 
mkdir -p $SLURM_SUBMIT_DIR ; 
# cp -r ./ $SLURM_SUBMIT_DIR
# cd $SLURM_SUBMIT_DIR || exit
##  use unique mcr cache directory on /scratch:
echo ramdisk folder created!

run='G_eq'
#outdir='/scratch/project/open-30-66/geq_data/'

## Loop through cores
for n in `seq 1 $SLURM_JOB_NUM_NODES`; do
 for i in `seq 1 36`; do
    node=$((n-1))
    # node=$(( $SLURM_NODEID  ))
    proc=$(( 36*(node)+i ))
	echo starting job $proc ... using node $node
	##  use unique mcr cache directory on /ramdisk:
	export MCR_CACHE_ROOT=`mktemp -d /ramdisk/$SLURM_JOB_ID/mcr.XXXXXXXXXX` 
	echo MCR_CACHE_FOLDER created at $MCR_CACHE_ROOT
	# cd $PBS_O_WORKDIR
	tasknumber=$(( i-1 ))
	srun -N 1 -n 1 --exclusive ./$run $proc 144 &
 done
done
 


 echo All done! ..
 wait
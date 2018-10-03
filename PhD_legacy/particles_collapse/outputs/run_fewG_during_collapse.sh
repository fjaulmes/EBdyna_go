#PBS -lwalltime=48:00:00
                         # 48 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores8:ppn=2
                         # 1 node for this job
module load mcr
# create unique mcr cache directory on /scratch:
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX` 
cd $HOME/MATLAB/JET_85383_ms/particles_collapse/outputs
./GT_few_evolution_collapse

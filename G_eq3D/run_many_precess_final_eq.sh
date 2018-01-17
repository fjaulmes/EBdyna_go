#PBS -lwalltime=28:00:00
                         # 28 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores8:ppn=2
                         # 1 node for this job
module load mcr
# create unique mcr cache directory on /scratch:
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX` 
cd $HOME/MATLAB/JET_060513/particles_equilibrium
./GT_many_precession_final_evolution_eq

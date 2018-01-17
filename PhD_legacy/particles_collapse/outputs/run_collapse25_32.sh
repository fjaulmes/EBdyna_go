#PBS -lwalltime=118:00:00
                         # 118 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=2
                         # 1 node for this job
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
module load mcr
( 
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 0 ./G_evolution_collapse 25
) &
(
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 1 ./G_evolution_collapse 26
) &
(
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 2 ./G_evolution_collapse 27
) &
(
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 3 ./G_evolution_collapse 28
) &
( 
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 4 ./G_evolution_collapse 29
) &
(
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 5 ./G_evolution_collapse 30
) &
(
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 6 ./G_evolution_collapse 31
) &
(
cd $HOME/MATLAB/AUG_31557_2p25_new/particles_collapse/outputs
taskset -c 7 ./G_evolution_collapse 32
) &
wait
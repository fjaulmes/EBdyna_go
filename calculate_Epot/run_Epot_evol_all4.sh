#PBS -lwalltime=22:00:00
                         # 22 hours wall-clock 
                         # time allowed for this job 
#PBS -lnodes=1
                         # 1 node for this job 
module load mcr 
# create unique mcr cache directory on /scratch:
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX` 
( 
cd $HOME/MATLAB/JET_180915h/calculate_Epot
taskset -c 0 ./calculate_E_potential_phi_values_evol13 
) &
(
cd $HOME/MATLAB/JET_180915h/calculate_Epot
taskset -c 1 ./calculate_E_potential_phi_values_evol14 
) &
(
cd $HOME/MATLAB/JET_180915h/calculate_Epot
taskset -c 2 ./calculate_E_potential_phi_values_evol15 
) &
(
cd $HOME/MATLAB/JET_180915h/calculate_Epot
taskset -c 3 ./calculate_E_potential_phi_values_evol16 
) &
wait



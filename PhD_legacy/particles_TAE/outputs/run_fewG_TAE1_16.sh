#PBS -lwalltime=120:00:00
                         # 120 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=1
                         # 1 node for this job
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
module load mcr
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 0 ./G_evolution_TAE  1
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 1 ./G_evolution_TAE  2
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 2 ./G_evolution_TAE  3
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 3 ./G_evolution_TAE  4
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 4 ./G_evolution_TAE  5
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 5 ./G_evolution_TAE  6
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 6 ./G_evolution_TAE  7
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 7 ./G_evolution_TAE  8
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 8 ./G_evolution_TAE  9
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 9 ./G_evolution_TAE  10
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 10 ./G_evolution_TAE  11
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 11 ./G_evolution_TAE  12
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 12 ./G_evolution_TAE  13
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 13 ./G_evolution_TAE  14
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 14 ./G_evolution_TAE  15
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 15 ./G_evolution_TAE  16
) &
wait



#PBS -lwalltime=120:00:00
                         # 120 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=1
                         # 1 node for this job
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
module load mcr
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 0 ./G_evolution_TAE  17
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 1 ./G_evolution_TAE  18
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 2 ./G_evolution_TAE  19
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 3 ./G_evolution_TAE  20
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 4 ./G_evolution_TAE  21
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 5 ./G_evolution_TAE  22
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 6 ./G_evolution_TAE  23
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 7 ./G_evolution_TAE  24
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 8 ./G_evolution_TAE  25
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 9 ./G_evolution_TAE  26
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 10 ./G_evolution_TAE  27
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 11 ./G_evolution_TAE  28
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 12 ./G_evolution_TAE  29
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 13 ./G_evolution_TAE  30
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 14 ./G_evolution_TAE  31
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 15 ./G_evolution_TAE  32
) &
wait



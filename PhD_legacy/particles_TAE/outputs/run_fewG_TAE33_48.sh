#PBS -lwalltime=120:00:00
                         # 120 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=1
                         # 1 node for this job
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
module load mcr
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 0 ./G_evolution_TAE  33
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 1 ./G_evolution_TAE  34
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 2 ./G_evolution_TAE  35
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 3 ./G_evolution_TAE  36
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 4 ./G_evolution_TAE  37
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 5 ./G_evolution_TAE  38
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 6 ./G_evolution_TAE  39
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 7 ./G_evolution_TAE  40
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 8 ./G_evolution_TAE  41
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 9 ./G_evolution_TAE  42
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 10 ./G_evolution_TAE  43
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 11 ./G_evolution_TAE  44
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 12 ./G_evolution_TAE  45
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 13 ./G_evolution_TAE  46
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 14 ./G_evolution_TAE  47
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 15 ./G_evolution_TAE  48
) &
wait



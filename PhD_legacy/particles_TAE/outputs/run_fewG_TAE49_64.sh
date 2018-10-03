#PBS -lwalltime=120:00:00
                         # 120 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=1
                         # 1 node for this job
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
module load mcr
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 0 ./G_evolution_TAE  49
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 1 ./G_evolution_TAE  50
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 2 ./G_evolution_TAE  51
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 3 ./G_evolution_TAE  52
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 4 ./G_evolution_TAE  53
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 5 ./G_evolution_TAE  54
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 6 ./G_evolution_TAE  55
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 7 ./G_evolution_TAE  56
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 8 ./G_evolution_TAE  57
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 9 ./G_evolution_TAE  58
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 10 ./G_evolution_TAE  59
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 11 ./G_evolution_TAE  60
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 12 ./G_evolution_TAE  61
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 13 ./G_evolution_TAE  62
) &
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 14 ./G_evolution_TAE  63
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/particles_TAE/outputs
taskset -c 15 ./G_evolution_TAE  64
) &
wait



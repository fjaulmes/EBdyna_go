#PBS -lwalltime=120:00:00
                         # 120 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=2
                         # 1 node for this job
module load mcr
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
( 
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 0 ./map_TAE_fields1
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 1 ./map_TAE_fields2
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 2 ./map_TAE_fields3
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 3 ./map_TAE_fields4
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 4 ./map_TAE_fields5
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 5 ./map_TAE_fields6
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 6 ./map_TAE_fields7
) &
(
cd $HOME/MATLAB/ITPA_benchmark_v3/build_TAE_maps
taskset -c 7 ./map_TAE_fields8
) &
wait



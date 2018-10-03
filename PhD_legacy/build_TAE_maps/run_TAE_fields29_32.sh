#PBS -lwalltime=113:00:00
                         # 113 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=2
                         # 1 node for this job
module load mcr
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX`
( 
cd $HOME/MATLAB/ITPA_benchmark/build_TAE_maps
taskset -c 0 ./map_TAE_fields29
) &
(
cd $HOME/MATLAB/ITPA_benchmark/build_TAE_maps
taskset -c 1 ./map_TAE_fields30
) &
(
cd $HOME/MATLAB/ITPA_benchmark/build_TAE_maps
taskset -c 2 ./map_TAE_fields31
) &
(
cd $HOME/MATLAB/ITPA_benchmark/build_TAE_maps
taskset -c 3 ./map_TAE_fields32
) &
wait



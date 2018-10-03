#PBS -lwalltime=89:00:00
                         # 89 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:cores16:ppn=1
                         # 1 node for this job
module load mcr
# create unique mcr cache directory on /scratch:
export MCR_CACHE_ROOT=`mktemp -d /scratch/mcr.XXXXXXXXXX` 
( 
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 0 ./GT_many_precession_post_evolution_eq1
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 1 ./GT_many_precession_post_evolution_eq2
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 2 ./GT_many_precession_post_evolution_eq3
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 3 ./GT_many_precession_post_evolution_eq4
) &
( 
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 4 ./GT_many_precession_post_evolution_eq5
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 5 ./GT_many_precession_post_evolution_eq6
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 6 ./GT_many_precession_post_evolution_eq7
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 7 ./GT_many_precession_post_evolution_eq8
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 8 ./GT_many_precession_post_evolution_eq9
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 9 ./GT_many_precession_post_evolution_eq10
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 10 ./GT_many_precession_post_evolution_eq11
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 11 ./GT_many_precession_post_evolution_eq12
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 12 ./GT_many_precession_post_evolution_eq13
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 13 ./GT_many_precession_post_evolution_eq14
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 14 ./GT_many_precession_post_evolution_eq15
) &
(
cd $HOME/MATLAB/AUG_30809/particles_equilibrium
taskset -c 15 ./GT_many_precession_post_evolution_eq16
) &
wait


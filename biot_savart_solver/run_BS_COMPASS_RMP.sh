#PBS -lwalltime=120:00:00
                         # 2 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:mem64gb
                         # 1 node for this job
						 #PBS -S /bin/bash
module load stopos
## FILL IN NAME OF SCRIPT
run='./BS_COMPASS_RMP_main'
# Number of parts in this pool
m=`stopos status -p pool_RMP 2>&1 | grep added | awk '{print $3;}'`
echo "Number of parts $m"

## Make directories for execution and output
cd $PBS_O_WORKDIR
#rm -r ./output

## Copy to TMPDIR
mkdir -p $TMPDIR/execution_folder/input
mkdir -p $TMPDIR/execution_folder/output
cp -r  ../data_tokamak $TMPDIR
cp -r  $PBS_O_WORKDIR/input $TMPDIR/execution_folder
cp $PBS_O_WORKDIR/$run $TMPDIR/execution_folder
cp $PBS_O_WORKDIR/task_script_RMP.sh $TMPDIR/execution_folder


## Nice information about nodes
# n = total number of tasks on a single node
n=`cat /proc/cpuinfo | grep processor | wc -l`
echo start of job in directory $PBS_O_WORKDIR
echo number of cores is $n

## Parallel execution
cd $TMPDIR/execution_folder
for i in `seq 1 $n`; do
   ./task_script_RMP.sh $run &
done
wait

## Email notification
echo "YEAH, your job $PBS_JOBID of file: $PBS_JOBNAME finished at `date`!" | mail $USER -s "FINISHED: Job $PBS_JOBID"
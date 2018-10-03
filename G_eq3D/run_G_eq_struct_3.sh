#PBS -lwalltime=120:00:00
                         # 60 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:mem64gb
                         # 1 node for this job
						 #PBS -S /bin/bash
module load stopos

## FILL IN NAME OF SCRIPT THAT NEEDS TO BE RUN BEFORE SIMULATION SCRIPT
prerun='./load_RAM'
## FILL IN NAME OF SCRIPT
run='./G_eq'

# Number of parts in this pool
m=`stopos status -p pool3 2>&1 | grep added | awk '{print $3;}'`
echo "Number of parts $m"

## Make directories for execution and output
cd $PBS_O_WORKDIR
mkdir -p $TMPDIR/execution_folder/input
mkdir -p $TMPDIR/execution_folder/output

# Remove any (previous) shared RAM information
rm -f ./busy.txt
rm -f ./shared_memory_keys.mat

## Copy to TMPDIR
cp -r  ../data_tokamak $TMPDIR
cp -r  $PBS_O_WORKDIR/input $TMPDIR/execution_folder
cp -r  $PBS_O_WORKDIR/output_14xx/* $TMPDIR/execution_folder/output/
cp $PBS_O_WORKDIR/$prerun $TMPDIR/execution_folder
cp $PBS_O_WORKDIR/$run $TMPDIR/execution_folder
cp $PBS_O_WORKDIR/task_script_3.sh $TMPDIR/execution_folder

## Nice information about nodes
# n = total number of tasks on a single node
n=`cat /proc/cpuinfo | grep processor | wc -l`
echo start of job in directory $PBS_O_WORKDIR
echo number of cores is $n

## load_RAM by TMPDIR
cd $TMPDIR/execution_folder
module load mcr/2016a
$prerun	

## Parallel execution
for i in `seq 1 $n`; do
   ./task_script_3.sh $run &
done
wait

$prerun 1	#Remove general RAM variables

## Email notification
echo "YEAH, your job $PBS_JOBID of file: $PBS_JOBNAME finished at `date`!" | mail $USER -s "FINISHED: Job $PBS_JOBID"
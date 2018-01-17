#PBS -lwalltime=120:00:00
                         # 120 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:mem64gb
                         # 1 node for this job
						 #PBS -S /bin/bash
module load stopos

## FILL IN NAME OF SCRIPT THAT NEEDS TO BE RUN BEFORE SIMULATION SCRIPT
prerun="./load_RAM_$ID"
## FILL IN NAME OF SCRIPT
run="./G_eq_$ID"

# Number of parts in this pool
m=`stopos status -p pool$ID 2>&1 | grep added | awk '{print $3;}'`
echo "Number of parts $m"

## Make directories for execution and output
cd $PBS_O_WORKDIR
mkdir -p $TMPDIR/execution_folder/input
mkdir -p $TMPDIR/execution_folder/output

# Remove any (previous) shared RAM information
rm -f ./busy.txt
rm -f ./shared_memory_keys.mat

## Copy to TMPDIR
cp $PBS_O_WORKDIR/output/*_prec.mat $TMPDIR/execution_folder/output/
cp -r  $PBS_O_WORKDIR/output_$ID/* $TMPDIR/execution_folder/output/
cp $PBS_O_WORKDIR/$prerun $TMPDIR/execution_folder
cp $PBS_O_WORKDIR/$run $TMPDIR/execution_folder
cp $PBS_O_WORKDIR/task_script_i.sh $TMPDIR/execution_folder
cp -r  ../data_tokamak $TMPDIR
cp -r  $PBS_O_WORKDIR/input/*.mat $TMPDIR/execution_folder/input/
#cp $PBS_O_WORKDIR/ba_interp2.mexa64 $TMPDIR/execution_folder
#cp $PBS_O_WORKDIR/ba_interp3.mexa64 $TMPDIR/execution_folder
#cp $PBS_O_WORKDIR/sharedmatrix.mexa64 $TMPDIR/execution_folder

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
   ./task_script_i.sh $run $ID &
done
wait

$prerun 1	#Remove general RAM variables

## notification
echo "PBS info: your job $PBS_JOBID of file: $PBS_JOBNAME finished at `date`!" 

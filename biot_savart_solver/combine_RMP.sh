#PBS -lwalltime=04:00:00
                         # 2 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:mem64gb
                         # 1 node for this job
						 #PBS -S /bin/bash
cd $PBS_O_WORKDIR

## FILL IN NAME OF SCRIPT
run='./BS_RMP_combine_output'

module load matlab/2017b
mcc -m $run.m
rm -f run_BS_RMP_combine_output.sh
rm -f *.txt
rm -f *.log

## Copy to TMPDIR
mkdir -p $TMPDIR/execution_folder/input
mkdir -p $TMPDIR/execution_folder/output
cp -r  ../data_tokamak $TMPDIR
cp -r  $PBS_O_WORKDIR/input $TMPDIR/execution_folder
cp -r  $PBS_O_WORKDIR/output_RMP $TMPDIR/execution_folder
cp -r  $PBS_O_WORKDIR/output/RMP*single_coil* $TMPDIR/execution_folder/output/
cp $run $TMPDIR/execution_folder

## Execute on scratch
cd $TMPDIR/execution_folder

## Nice information about nodes
# n = total number of tasks on a single node
module load mcr/v90
$run 

## Copy back output
cp -r $TMPDIR/execution_folder/output $PBS_O_WORKDIR

## Email notification
echo "YEAH, your job $PBS_JOBID of file: $PBS_JOBNAME finished at `date`!" | mail $USER -s "FINISHED: Job $PBS_JOBID"
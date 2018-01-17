#PBS -lwalltime=10:00:00
                         # 60 hours wall-clock 
                         # time allowed for this job
#PBS -lnodes=1:mem64gb
                         # 1 node for this job
						 #PBS -S /bin/bash
## FILL IN NAME OF SCRIPT
run='./combine_output_eq'

## Make directories for execution and output
cd $PBS_O_WORKDIR
mkdir -p $PBS_O_WORKDIR/output

module load matlab/2016a
cp  ../01\ general_functions/* ./
mcc -m combine_output_eq.m
rm -f run_combine_output_eq.sh
rm -f *.txt
rm -f *.log

## Copy to TMPDIR
mkdir -p $TMPDIR/execution_folder/output
mkdir -p $TMPDIR/execution_folder/input

cp -r  ../data_tokamak $TMPDIR
cp $PBS_O_WORKDIR/$run $TMPDIR/execution_folder

## export ID=19xx
echo the ID of simulation folder where files will be combined is $ID


cp -r  $PBS_O_WORKDIR/output_$ID $TMPDIR/execution_folder

module load mcr/2016a
cd $TMPDIR/execution_folder
$run $ID prec



## Copy back output
cp -r $TMPDIR/execution_folder/input $PBS_O_WORKDIR
cp -r $TMPDIR/execution_folder/output $PBS_O_WORKDIR

## Email notification
echo "YEAH, your job $PBS_JOBID of file: $PBS_JOBNAME finished at `date`!" 

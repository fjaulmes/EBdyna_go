#!/bin/bash
module load stopos
module load mcr/2016a
if [ "x$1" == "x" ]; then
   echo "Please use binary name on command line!"
   exit 1
fi
run=$1
m=`stopos status -p pool_RMP_15611 2>&1 | grep added | awk '{print $3;}'`
STOPOS_RC="OK"
while [ "x$STOPOS_RC" == "xOK" ]; do
	stopos next -m -p pool_RMP_15611
	if [ $STOPOS_COMMITTED -gt 0 ]; then
		STOPOS_RC="DONE"
	else
		if [ -f $PBS_O_WORKDIR/output_RMP/*_process$STOPOS_VALUE.mat ]; then
			stopos remove -p pool_RMP_15611
		else
			$run $STOPOS_VALUE $m
			stopos remove -p pool_RMP_15611
			mv $TMPDIR/execution_folder/output/*_process$STOPOS_VALUE.mat* $PBS_O_WORKDIR/output_RMP
		fi
	fi
done

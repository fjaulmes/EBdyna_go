#!/bin/bash
module load stopos
module load mcr/2016a
if [ "x$1" == "x" ]; then
   echo "Please use binary name on command line!"
   exit 1
fi
run=$1
m=`stopos status -p pool3 2>&1 | grep added | awk '{print $3;}'`
STOPOS_RC="OK"
while [ "x$STOPOS_RC" == "xOK" ]; do
	stopos next -m -p pool3
    if [ $STOPOS_COMMITTED -gt 0 ]; then
		STOPOS_RC="DONE"
	else
		if [ -f $PBS_O_WORKDIR/output_3/*_full_process$STOPOS_VALUE.mat* ]; then
			stopos remove -p pool3
		else
			$run $STOPOS_VALUE $m
			stopos remove -p pool3
			mv $TMPDIR/execution_folder/output/*_process$STOPOS_VALUE.mat* $PBS_O_WORKDIR/output_3
		fi
	fi
done

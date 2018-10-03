#!/bin/bash
module load stopos
module load mcr/2016a
run=$1
ID=$2
m=`stopos status -p pool$ID 2>&1 | grep added | awk '{print $3;}'`
STOPOS_RC="OK"
while [ "x$STOPOS_RC" == "xOK" ]; do
	stopos next -m -p pool$ID
    if [ $STOPOS_COMMITTED -gt 0 ]; then
		STOPOS_RC="DONE"
	else
		if [ -f $PBS_O_WORKDIR/output_$ID/*_full_process$STOPOS_VALUE.mat* ]; then
			stopos remove -p pool$ID
		else
			$run $STOPOS_VALUE $m
			stopos remove -p pool$ID
			rm -f $TMPDIR/execution_folder/output/G_eq_*xx_prec_process$STOPOS_VALUE.mat
			mv $TMPDIR/execution_folder/output/*_process$STOPOS_VALUE.mat* $PBS_O_WORKDIR/output_$ID
		fi
	fi
done

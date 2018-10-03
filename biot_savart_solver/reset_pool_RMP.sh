#!/bin/bash
#
m=192
module load stopos
stopos purge -p pool_RMP_15611
stopos create -p pool_RMP_15611
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool_RMP_15611 add parmset.$$
rm -f parmset.$$

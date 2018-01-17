#!/bin/bash
#
m=1000
module load stopos
stopos purge -p pool_TFR
stopos create -p pool_TFR
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool_TFR add parmset.$$
rm -f parmset.$$

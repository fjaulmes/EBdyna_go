#!/bin/bash
#
m=100
module load stopos
stopos purge -p pool_RMP
stopos create -p pool_RMP
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool_RMP add parmset.$$
rm -f parmset.$$

#!/bin/bash
#
m=1000
module load stopos
stopos purge -p pool1
stopos create -p pool1
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool1 add parmset.$$
rm -f parmset.$$

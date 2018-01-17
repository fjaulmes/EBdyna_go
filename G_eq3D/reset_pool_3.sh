#!/bin/bash
#
m=300
module load stopos
stopos purge -p pool3
stopos create -p pool3
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool3 add parmset.$$
rm -f parmset.$$

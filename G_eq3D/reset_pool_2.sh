#!/bin/bash
#
m=300
module load stopos
stopos purge -p pool2
stopos create -p pool2
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool2 add parmset.$$
rm -f parmset.$$

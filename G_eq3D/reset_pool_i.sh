#!/bin/bash
#
m=200
module load stopos
stopos purge -p pool$ID
stopos create -p pool$ID
rm -f parmset.$$
for i in `seq 1 $m`; do
   echo $i >> parmset.$$
done
stopos -p pool$ID add parmset.$$
rm -f parmset.$$

#!/bin/bash

proc=$((SLURM_PROCID + 1))
total=$SLURM_NTASKS

# Ensure parent folder exists
SLURM_SUBMIT_DIR=/ramdisk/$SLURM_JOB_ID
mkdir -p $SLURM_SUBMIT_DIR

# Unique MCR cache per task
export MCR_CACHE_ROOT=$SLURM_SUBMIT_DIR/mcr_$proc
mkdir -p $MCR_CACHE_ROOT

echo "Running Kk-task $proc of $total with MCR cache at $MCR_CACHE_ROOT"

# Run the MATLAB compiled executable
./G_eq_203301_30L_1DF06 $proc $total
#!/bin/sh

DIR=$1

# Change to the example directory
cd $DIR

# Set the vtune config script.
source /projects/partisn/vtune_amplifier_xe/amplxe-vars.sh

# Set the thread policy
export KMP_AFFINITY=SCATTER

# Set the number of threads
export OMP_NUM_THREADS=240

# Set the location of the intel shared libraries
export LD_LIBRARY_PATH=/projects/partisn/composerxe/lib/mic 

# Unlimit the stack size
ulimit -s unlimited

# Launch the executable
../../build-mic-15/src/hc/pececillo-HC ex1.json


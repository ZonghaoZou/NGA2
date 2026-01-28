#!/bin/bash

for i in {0..9}; do
    # Directly edit the input file for this run
    sed "s/^Meshratio.*/Meshratio : $i/" input > input.tmp
    mv input.tmp input

    # Create a folder for this case
    mkdir -p case_$i
    cp input case_$i/

    # Run the solver inside the case folder
    echo "Running case $i with Meshratio=$i"
    cd case_$i
    mpiexec -n 6 ../nga.dp.gnu.opt.mpi.exe -i input -v 2 > run.log
    cd ..
done
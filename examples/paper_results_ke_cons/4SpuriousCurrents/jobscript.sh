#!/bin/bash

NUM_RUNS=30

for((i=1;i<=NUM_RUNS;i++))
do 
	echo "Running simulation $i..."
	mpiexec -n 1 ./nga.dp.gnu.opt.mpi.exe -i input -v 2
done

echo "All simulations completed."


#!/bin/bash

# --- Configuration Arrays ---
# Aligned indices (0 to 3) for Mesh, Partition, and Cores
NX_VALUES=(64 32 16 8)
PARTITION_VALUES=("6 1 1" "4 1 1" "2 1 1" "1 1 1")
CORE_VALUES=(6 4 2 1)

# Timesteps for each Laplace number (Aligned with NX_VALUES above)
# Index 0 (nx=64) gets smallest dt, Index 3 (nx=8) gets biggest dt
DT_LA_12e2=("1.2e-3" "3.4e-3" "9.7e-3" "2.7e-2")
DT_LA_12e4=("1.2e-2" "3.4e-2" "9.7e-2" "2.7e-1")
DT_LA_12e6=("1.2e-1" "3.4e-1" "9.7e-1" "2.7")

LA_VALUES=("1.2e+2" "1.2e+4" "1.2e+6")

# Number of times to repeat the simulation for each case
NUM_RUNS=30

# Backup original input
cp input input.bak

# --- Main Parameter Loops ---
for la in "${LA_VALUES[@]}"; do
    
    # Select the correct DT array based on Laplace number
    if [ "$la" == "1.2e+2" ]; then
        current_dt_list=("${DT_LA_12e2[@]}")
    elif [ "$la" == "1.2e+4" ]; then
        current_dt_list=("${DT_LA_12e4[@]}")
    elif [ "$la" == "1.2e+6" ]; then
        current_dt_list=("${DT_LA_12e6[@]}")
    fi

    # Loop through mesh sizes
    for i in "${!NX_VALUES[@]}"; do
        # Extract parameters
        nx=${NX_VALUES[$i]}
        partition=${PARTITION_VALUES[$i]}
        cores=${CORE_VALUES[$i]}
        dt=${current_dt_list[$i]}

        # Define run directory
        run_dir="Run_La${la}_nx${nx}"
        
        echo "======================================================"
        echo "Setting up Case: La=$la | nx=$nx | Part=($partition) | Cores=$cores | dt=$dt"
        echo "Directory: $run_dir"

        # Prepare directory
        mkdir -p "$run_dir"
        cp nga.dp.gnu.opt.mpi.exe "$run_dir/"
        cp input "$run_dir/"
        
        # Enter the directory
        cd "$run_dir" || exit

        # --- Modify input file (Only needs to be done once per folder) ---
        sed -i "s/Partition : .*/Partition : $partition/" input
        sed -i "s/nx : .*/nx : $nx/" input
        sed -i "s/Laplace number:.*/Laplace number:               $la/" input
        sed -i "s/Max timestep size : .*/Max timestep size : $dt/" input

        # --- 30 Runs Loop ---
        echo "  -> Starting batch of $NUM_RUNS runs..."
        
        for ((k=1; k<=NUM_RUNS; k++)); do
            echo "     Running iteration $k of $NUM_RUNS..."
            
            # Using append (>>) for log so we keep history of all 30 runs
            echo "----- Run $k -----" >> log.txt
            mpiexec -n "$cores" ./nga.dp.gnu.opt.mpi.exe -i input -v 2 >> log.txt 2>&1
        done
        
        echo "  -> Batch finished."
        
        # Return to main folder
        cd ..
    done
done

echo "======================================================"
echo "All study simulations completed."
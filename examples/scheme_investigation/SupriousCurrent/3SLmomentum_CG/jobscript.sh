#!/bin/bash

# --- Configuration Arrays (Reversed: Low Nodes -> High Nodes) ---
# Index 0: nx=8,  Cores=1, Part=1 1 1
# Index 3: nx=64, Cores=6, Part=6 1 1

NX_VALUES=(8 16 32 64)
PARTITION_VALUES=("1 1 1" "2 1 1" "4 1 1" "6 1 1")
CORE_VALUES=(1 2 4 6)

# Timesteps (Reversed to match NX_VALUES)
# Index 0 (nx=8) gets Biggest DT
# Index 3 (nx=64) gets Smallest DT
DT_LA_12e2=("2.7e-2" "9.7e-3" "3.4e-3" "1.2e-3")
DT_LA_12e4=("2.7e-1" "9.7e-2" "3.4e-2" "1.2e-2")
DT_LA_12e6=("2.7"    "9.7e-1" "3.4e-1" "1.2e-1")

LA_VALUES=("1.2e+2" "1.2e+4" "1.2e+6")

NUM_RUNS=30

# Backup original input
cp input input.bak

# --- Main Loops ---
for la in "${LA_VALUES[@]}"; do
    
    # Select DT list based on Laplace number
    if [ "$la" == "1.2e+2" ]; then
        current_dt_list=("${DT_LA_12e2[@]}")
    elif [ "$la" == "1.2e+4" ]; then
        current_dt_list=("${DT_LA_12e4[@]}")
    elif [ "$la" == "1.2e+6" ]; then
        current_dt_list=("${DT_LA_12e6[@]}")
    fi

    for i in "${!NX_VALUES[@]}"; do
        nx=${NX_VALUES[$i]}
        partition=${PARTITION_VALUES[$i]}
        cores=${CORE_VALUES[$i]}
        dt=${current_dt_list[$i]}

        run_dir="Run_La${la}_nx${nx}"
        
        echo "======================================================"
        echo "Setting up: La=$la | nx=$nx | Part=($partition) | Cores=$cores | dt=$dt"

        mkdir -p "$run_dir"
        cp nga.dp.gnu.opt.mpi.exe "$run_dir/"
        cp input "$run_dir/"
        
        cd "$run_dir" || exit

        # --- FIX: Flexible SED Commands for Mac ---
        # matches line starting with "Partition" followed by anything until ":"
        sed -i '' "s/^Partition.*:.*/Partition : $partition/" input
        sed -i '' "s/^nx.*:.*/nx : $nx/" input
        sed -i '' "s/^Laplace number.*:.*/Laplace number:               $la/" input
        sed -i '' "s/^Max timestep size.*:.*/Max timestep size : $dt/" input

        # --- Verify the change (Optional check printed to screen) ---
        echo "   [Check] Input file partition line is now:"
        grep "Partition" input

        # --- 30 Runs Loop ---
        echo "  -> Starting batch of $NUM_RUNS runs..."
        
        for ((k=1; k<=NUM_RUNS; k++)); do
            # Appending to log so we don't overwrite previous runs
            echo "----- Run $k -----" >> log.txt
            mpiexec -n "$cores" ./nga.dp.gnu.opt.mpi.exe -i input -v 2 >> log.txt 2>&1
        done
        
        cd ..
    done
done

echo "All simulations completed."
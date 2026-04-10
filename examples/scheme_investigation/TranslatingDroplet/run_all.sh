#!/bin/bash
# =========================================================================
# run_all.sh — Run translating droplet test across all mesh resolutions
#              and all three NGA2 schemes
#
# D/Δ = 6.4, 12.8, 25.6, 51.2, 102.4
# nx  = 16,  32,   64,   128,   256
#
# Following the spurious current jobscript.sh pattern:
#   1. Copy base input file
#   2. Use sed to modify resolution-specific params
#   3. Run in subdirectory
# =========================================================================
set -e

BASEDIR="$(cd "$(dirname "$0")" && pwd)"
SCHEMES=("1DefaultNGA2_SG" "2SpatialTemporalKEcons_SG" "3SLmomentum_CG")

# Resolution parameters
# NX_VALUES=(16 32 64 128 256)
# PARTITION_VALUES=("2 2 1" "2 2 1" "2 2 1" "3 2 1" "3 2 1")
# CORE_VALUES=(4 4 4 6 6)
# DT_VALUES=("1e-5" "5e-6" "1e-6" "3e-7" "1e-7")

NX_VALUES=(16 32 64)
PARTITION_VALUES=("2 2 1" "2 2 1" "3 2 1")
CORE_VALUES=(4 4 6)
DT_VALUES=("1.8e-6" "6.5e-7" "2.3e-7")

# =========================================================================
# Phase 1: Compile all three schemes
# =========================================================================
echo "============================================"
echo "  Phase 1: Compiling all schemes"
echo "============================================"
for scheme in "${SCHEMES[@]}"; do
    echo "--- Compiling ${scheme} ---"
    cd "${BASEDIR}/${scheme}"
    make -j6 2>&1 | tail -3
    echo "    Done."
done

# =========================================================================
# Phase 2: Run all resolutions for all schemes
# =========================================================================
echo ""
echo "============================================"
echo "  Phase 2: Running simulations"
echo "============================================"

for i in "${!NX_VALUES[@]}"; do
    nx=${NX_VALUES[$i]}
    partition=${PARTITION_VALUES[$i]}
    cores=${CORE_VALUES[$i]}
    dt=${DT_VALUES[$i]}
    
    echo ""
    echo "======== Resolution: nx=${nx} (D/Delta=$(echo "scale=1; ${nx}/2.5" | bc)) ========"
    
    for scheme in "${SCHEMES[@]}"; do
        echo "  --- ${scheme}, nx=${nx} ---"
        cd "${BASEDIR}/${scheme}"
        
        run_dir="Run_nx${nx}"
        
        mkdir -p "$run_dir"
        cp nga.dp.gnu.opt.mpi.exe "$run_dir/"
        cp "${BASEDIR}/input" "$run_dir/input"
        
        cd "$run_dir" || exit
        
        # Modify input file with sed (Mac-compatible)
        sed -i '' "s/^Partition.*:.*/Partition : $partition/" input
        sed -i '' "s/^nx.*:.*/nx : $nx/" input
        sed -i '' "s/^ny.*:.*/ny : $nx/" input
        sed -i '' "s/^nz.*:.*/nz : $nx/" input
        sed -i '' "s/^Max timestep size.*:.*/Max timestep size : $dt/" input
        
        echo "      [Check] Partition: $(grep 'Partition' input)"
        echo "      [Check] nx: $(grep '^nx' input)"
        echo "      [Check] dt: $(grep 'Max timestep size' input)"
        
        echo "      mpiexec -n ${cores} ./nga.dp.gnu.opt.mpi.exe -i input -v 2"
        mpiexec -n ${cores} ./nga.dp.gnu.opt.mpi.exe -i input -v 2 >> log.txt 2>&1
        
        # Copy CSV result to consolidated location
        mkdir -p "${BASEDIR}/result/${scheme}"
        if [ -f "${nx}.csv" ]; then
            cp "${nx}.csv" "${BASEDIR}/result/${scheme}/"
            echo "      -> Saved ${nx}.csv"
        fi
        
        cd "${BASEDIR}/${scheme}"
        echo "      Done."
    done
done

echo ""
echo "============================================"
echo "  All simulations complete!"
echo "  Results in: ${BASEDIR}/result/"
echo "============================================"

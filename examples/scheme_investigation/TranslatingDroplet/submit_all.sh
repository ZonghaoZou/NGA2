#!/bin/bash
# =========================================================================
# submit_all.sh — Set up and submit translating droplet SLURM jobs
#
# Directory structure:
#   TranslatingDrop/
#   ├── input              (base input file)
#   ├── submit_all.sh      (this script)
#   ├── 1/{16,32,64,128,256}/   (Default NGA2)
#   ├── 2/{16,32,64,128,256}/   (KE Conservative)
#   └── 3/{16,32,64,128,256}/   (SL Momentum)
#
# Each resolution directory gets: source, GNUmakefile, input, job.script
# Then compiles and submits.
# =========================================================================

base_dir=$(pwd)

# --- Configuration ---
schemes=("1" "2" "3")
resolutions=(16 32 64 128 256)

#                nx=16       nx=32       nx=64       nx=128      nx=256
PART_VALUES=(   "3 2 2"     "3 2 2"     "4 3 2"     "4 3 2"     "6 5 4"  )
NTASKS=(        12          12          24          24          120      )
NODES=(         1           1           2           2           10       )
NTASKS_NODE=(   12          12          12          12          12       )
DT_VALUES=(     "1.8e-6"    "6.0e-7"    "2.0e-7"    "6.7e-8"    "2.2e-8" )
WALLTIME=(      "04:00:00"  "08:00:00"  "24:00:00"  "72:00:00"  "120:00:00")

SCHEME_NAMES=("DefaultNGA2" "KEcons" "SLmom")

# --- Step 1: Set up all directories ---
echo "============================================"
echo "  Setting up directories and input files"
echo "============================================"

for s_idx in "${!schemes[@]}"; do
    scheme=${schemes[$s_idx]}
    sname=${SCHEME_NAMES[$s_idx]}
    
    for r in "${!resolutions[@]}"; do
        nx=${resolutions[$r]}
        partition=${PART_VALUES[$r]}
        ntasks=${NTASKS[$r]}
        nodes=${NODES[$r]}
        ntasks_node=${NTASKS_NODE[$r]}
        dt=${DT_VALUES[$r]}
        wt=${WALLTIME[$r]}
        
        target_dir="${scheme}/${nx}"
        echo "  Setting up: Scheme ${scheme} | nx=${nx} -> ${target_dir}"
        
        # Create directory and copy build files + source from scheme base
        mkdir -p "${target_dir}/src"
        cp "${scheme}/GNUmakefile"   "${target_dir}/" 2>/dev/null
        cp "${scheme}/Make.package"  "${target_dir}/" 2>/dev/null
        cp "${scheme}/src/"*         "${target_dir}/src/" 2>/dev/null
        
        # Copy and modify input file
        cp input "${target_dir}/input"
        cd "${target_dir}"
        sed -i "s/^Partition.*:.*/Partition : $partition/" input
        sed -i "s/^nx.*:.*/nx : $nx/" input
        sed -i "s/^ny.*:.*/ny : $nx/" input
        sed -i "s/^nz.*:.*/nz : $nx/" input
        sed -i "s/^Max timestep size.*:.*/Max timestep size : $dt/" input
        
        # Create job.script
        cat > job.script <<JOBEOF
#!/bin/bash
#SBATCH -p normal
#SBATCH -t ${wt}
#SBATCH --nodes=${nodes}
#SBATCH --ntasks-per-node=${ntasks_node}
#SBATCH --exclude=c-3-[13,15]
#SBATCH --exclusive
#SBATCH -J TD_${sname}_${nx}
#SBATCH -o %j.out
#SBATCH -e %j.err
#SBATCH --mail-type=ALL
#SBATCH --mail-user=zz474@cornell.edu
mpiexec --mca io ^ompio nga.dp.gnu.opt.mpi.exe -i input -v 2
JOBEOF
        
        cd "$base_dir"
    done
done

echo ""
echo "============================================"
echo "  Compiling and submitting all jobs"
echo "============================================"

# --- Step 2: Compile and submit (following user's pattern) ---
for scheme in "${schemes[@]}"; do
    for r in "${!resolutions[@]}"; do
        nx=${resolutions[$r]}
        target_dir="${scheme}/${nx}"
        
        echo "========================================"
        echo "Processing: Scheme $scheme | Resolution $nx"
        
        if [ -d "$target_dir" ]; then
            cd "$target_dir" || exit
            
            echo "  > Compiling (make -j12)..."
            make -j12
            
            if [ $? -eq 0 ]; then
                echo "  > Compilation successful."
                
                if [ -f "job.script" ]; then
                    echo "  > Submitting to Slurm..."
                    sbatch job.script
                else
                    echo "  > ERROR: 'job.script' not found in $target_dir"
                fi
            else
                echo "  > ERROR: Make failed. Skipping job submission."
            fi
            
            cd "$base_dir" || exit
        else
            echo "  > WARNING: Directory '$target_dir' does not exist. Skipping."
        fi
    done
done

echo "========================================"
echo "All tasks completed."

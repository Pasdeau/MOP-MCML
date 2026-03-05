#!/bin/bash
#SBATCH --job-name=mcml_3d_gpu
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1
#SBATCH --mail-user=wenzheng.wang@lip6.fr
#SBATCH --mail-type=ALL
#SBATCH --time=01:30:00

# ----------------------------------------
# Initialize module system and load GCC 11 + CUDA
source /etc/profile.d/modules.sh
module purge
module load gcc/11 cuda/11.8
# ----------------------------------------

echo "============================================"
echo " MCML 3D GPU Version (Ronly, Tonly, RT)"
echo " $(date)"
echo "============================================"

# Common compilation flags
NVCC_FLAGS="-O3 -use_fast_math -arch=sm_80"
MCI_A="test_3d.mci"

compile_and_run() {
    local MODE=$1
    local LABEL=$2
    local OUTFILE="3d_gpu_${LABEL}_A.mco"
    
    echo "----------------------------------------"
    echo " [${LABEL}] Compiling (PD_MODE=${MODE})..."
    nvcc $NVCC_FLAGS -DPD_MODE=${MODE} -o mcml_3d_gpu_${LABEL} mcml_gpu.cu mcmlio_gpu.c
    
    if [ $? -ne 0 ]; then
        echo "Compilation failed for ${LABEL}!"
        return
    fi
    
    echo " [${LABEL}] Running simulation..."
    # The C code uses the output filename defined in the .mci file
    # We create a temporary .mci to override the output filename to match OUTFILE
    sed "s/test_3d\.mco/${OUTFILE}/g" ${MCI_A} > tmp_${LABEL}.mci
    
    ./mcml_3d_gpu_${LABEL} tmp_${LABEL}.mci
    rm -f tmp_${LABEL}.mci
    
    if [ -f "${OUTFILE}" ]; then
        echo " -> Output saved to ${OUTFILE}"
    else
        echo " -> ERROR: ${OUTFILE} was NOT produced."
    fi
}

compile_and_run 1 "Ronly"
compile_and_run 2 "Tonly"
compile_and_run 3 "RT"

echo ""
echo "==== Job Completed ===="
ls -lh *.mco

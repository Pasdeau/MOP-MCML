#!/bin/bash
#SBATCH --job-name=mcml_run
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1
#SBATCH --mail-user=wenzheng.wang@lip6.fr
#SBATCH --mail-type=ALL
#SBATCH --time=00:10:00

# Explicitly set library path since 'module' command might fail on compute nodes
# This path matches the login node's environment for cuda/11.8
export LD_LIBRARY_PATH=/usr/local/cuda-11.8/lib64:$LD_LIBRARY_PATH

cd ~/mcml_gpu_test

echo "Starting execution on $(hostname) at $(date)"

# Ensure binary exists
if [ ! -f ./mcml_gpu ]; then
    echo "Error: ./mcml_gpu not found!"
    exit 1
fi

# Clean up previous runs
rm -f *.mco

# Run simulations
echo "Running 3mm..."
./mcml_gpu 3mm_1500.mci

echo "Running 4mm..."
./mcml_gpu 4mm_1500.mci

echo "Running 5mm..."
./mcml_gpu 5mm_1500.mci

echo "All finished at $(date)"
ls -lh summary_*.csv

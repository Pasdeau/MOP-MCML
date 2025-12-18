#!/bin/bash
#SBATCH --job-name=mcml_batch
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1
#SBATCH --mail-user=wenzheng.wang@lip6.fr
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00

# Load environment
module purge
module load gcc/11 cuda/11.8

# Compile
cd ~/mcml_gpu_test
make

# Run simulations
echo "Starting batch simulation on $(hostname) at $(date)"

echo "Running 3mm_1500.mci..."
./mcml_gpu 3mm_1500.mci

echo "Running 4mm_1500.mci..."
./mcml_gpu 4mm_1500.mci

echo "Running 5mm_1500.mci..."
./mcml_gpu 5mm_1500.mci

echo "Batch simulation finished at $(date)"

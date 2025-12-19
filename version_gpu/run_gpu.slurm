#!/bin/bash
#SBATCH --job-name=mcml_gpu
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1
#SBATCH --mail-user=wenzheng.wang@lip6.fr
#SBATCH --mail-type=ALL
#SBATCH --time=02:00:00

# Load CUDA environment
module purge
module load gcc/11 cuda/11.8

# Run simulation
# Each .mci file will generate its own summary_<name>.csv
echo "Starting GPU simulation at $(date)"
./mcml_gpu input.mci
echo "Simulation finished at $(date)"

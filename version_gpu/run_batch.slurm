#!/bin/bash
#SBATCH --job-name=mcml_auto
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

# Recompile with latest code
cd ~/mcml_gpu_test
make clean
make

echo "Starting simulations on $(hostname) at $(date)"
echo "Using updated GPU code with auto-naming CSV feature"

# Clean up old .mco files
rm -f *.mco

# Run all 3 simulations - each will auto-generate its own CSV
./mcml_gpu 3mm_1500.mci
./mcml_gpu 4mm_1500.mci
./mcml_gpu 5mm_1500.mci

echo "All simulations finished at $(date)"
echo "Generated CSV files:"
ls -lh summary_*.csv

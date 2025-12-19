#!/bin/bash
#SBATCH --job-name=mcml_separate
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
make clean
make

# Run simulations separately and rename summary.csv after each
echo "Starting simulations on $(hostname) at $(date)"

# Clean up old .mco files and csv
rm -f *.mco summary*.csv

echo "Running 3mm_1500.mci..."
./mcml_gpu 3mm_1500.mci
if [ -f summary.csv ]; then
    mv summary.csv summary_3mm.csv
    echo "Generated summary_3mm.csv"
fi

echo "Running 4mm_1500.mci..."
./mcml_gpu 4mm_1500.mci
if [ -f summary.csv ]; then
    mv summary.csv summary_4mm.csv
    echo "Generated summary_4mm.csv"
fi

echo "Running 5mm_1500.mci..."
./mcml_gpu 5mm_1500.mci
if [ -f summary.csv ]; then
    mv summary.csv summary_5mm.csv
    echo "Generated summary_5mm.csv"
fi

echo "All simulations finished at $(date)"
echo "Generated CSV files:"
ls -lh summary_*.csv

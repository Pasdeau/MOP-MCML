#!/bin/bash
#SBATCH --job-name=mcml_gpu_cmp
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1
#SBATCH --time=01:30:00

module load gcc/11
module load cuda/11.8

echo "Compiling and running GPU R-only..."
nvcc -O3 mcml_gpu.cu mcmlio_gpu.c -o mcml_gpu_Ronly -DPD_MODE=1
./mcml_gpu_Ronly test_validate.mci
mv test_gpu.mco test_gpu_Ronly.mco

echo "Compiling and running GPU T-only..."
nvcc -O3 mcml_gpu.cu mcmlio_gpu.c -o mcml_gpu_Tonly -DPD_MODE=2
./mcml_gpu_Tonly test_validate.mci
mv test_gpu.mco test_gpu_Tonly.mco

echo "Compiling and running GPU R+T..."
nvcc -O3 mcml_gpu.cu mcmlio_gpu.c -o mcml_gpu_RT -DPD_MODE=3
./mcml_gpu_RT test_validate.mci
mv test_gpu.mco test_gpu_RT.mco

echo "GPU calculations finished."

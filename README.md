# MOP-MCML Project: Monte Carlo Multi-Layered Simulation

This project implements the Monte Carlo Multi-Layered (MCML) method for simulating light transport in multi-layered turbid media.
It relies on the method of **Mean Optical Path (MOP)** and includes both a standard **CPU version** and a **High-Performance GPU version**.

**Authors:**
*   Songlin Li – First modifications
*   Wenzheng Wang – Improvements and GPU port

**Supervisor:**
*   Prof. Sylvain Feruglio (Sorbonne Université, LIP6)

---

## 📂 Project Structure

*   `mcmlmain.c`, `mcmlgo.c`, `mcmlio.c`, etc.: **Standard CPU Source Code** (C Language).
*   `version_gpu/`: **GPU Source Code** (CUDA C/C++), including `mcml_gpu.cu`.
*   `.mci` files: Input configuration files defining layers and simulation parameters.
*   `.mco` files: Output result files.
*   `get_mop.m`: MATLAB script for visualization and analysis of results.

---

## 🚀 1. CPU Version (Standard)

This is the standard C implementation, suitable for small-scale simulations or debugging on systems without NVIDIA GPUs.

### Compilation (Windows / Visual Studio)
1.  **Open Project**: Locate the `.sln` file and open it with **Visual Studio 2022**.
    *   *Note*: VS 2019 projects may require "Retarget Projects" (right-click Solution -> Retarget) to upgrade the SDK version.
2.  **Build**: Select "Build Solution" (Ctrl+Shift+B). An `.exe` will be generated.

### Compilation (Linux / MacOS Terminal)
```bash
# Compile source files linked together
cl mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c
```
*(Or use `gcc` on Linux: `gcc -O2 -o mcml mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -lm`)*

---

## ⚡ 2. GPU Version (High Performance) & Comparison

The GPU version is designed to replace the CPU version for high-throughput simulations.

### 🆚 Performance Comparison

| Feature | CPU Version | GPU Version (CUDA) |
| :--- | :--- | :--- |
| **Simulation Time** | ~70.0 sec / run | **~0.4 sec / run** (on A100) |
| **Speedup** | 1x | **~175x Faster** |
| **Photon Transport** | Serial (One by one) | **Parallel** (Thousands at once) |
| **Precision** | Double Precision | Double Precision (Atomic Operations) |
| **Hardware** | Standard CPU | NVIDIA GPU (Compute Cap 5.0+) |

*> Note: Benchmark based on 2,000,000 photons per run.*

### GPU Version Location
Located in the `version_gpu/` folder. This version uses **CUDA** to parallelize photon transport.

### Prerequisites
*   **Hardware**: NVIDIA GPU (Compute Capability 5.0+, e.g., A100, V100, RTX A-series).
*   **Software**:
    *   CUDA Toolkit (e.g., v11.8).
    *   GCC (compatible version).

### Compilation
Navigate to the GPU directory and use `make`:
```bash
cd version_gpu
make
```
This generates the `mcml_gpu` executable.

### Usage
### Running a Simulation
```bash
./mcml_gpu <input_file.mci>
# Example:
./mcml_gpu input.mci
```

**Output**: Each `.mci` file automatically generates a corresponding CSV file:
- `3mm_1500.mci` → `summary_3mm_1500.csv`
- `4mm_1500.mci` → `summary_4mm_1500.csv`

### HPC / SLURM Submission
For cluster environments (e.g., LIP6 Convergence), use the provided SLURM scripts in `version_gpu/`:

**Single File Simulation:**
```bash
# Edit run_gpu.slurm to specify your input file
sbatch version_gpu/run_gpu.slurm
```

**Batch Processing (Multiple Files):**
```bash
# Runs multiple .mci files sequentially
sbatch version_gpu/run_batch.slurm
```

**SLURM Script Template:**
```bash
#!/bin/bash
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1
#SBATCH --mail-type=ALL

module purge
module load gcc/11 cuda/11.8

./mcml_gpu your_file.mci
```
*(The `run_batch.slurm` script can compile and run multiple simulations sequentially).*

---

## 📝 Input & Configuration (.mci Files)

The input file defines the simulation parameters.

**1. Photodiode (PD) Parameters** (2nd to last line):
Format (cm):
`Rx Ry Rl Tx Ty Tl`
*   **(Rx, Ry)**: Center position of Reflectance PD.
*   **Rl**: Side length of Reflectance PD.
*   **(Tx, Ty)**: Center position of Transmittance PD.
*   **Tl**: Side length of Transmittance PD.

**2. Light Source** (Last line):
Format:
`Type x y Param`
*   **Type**: `1` = Point, `2` = Gaussian, `3` = Flat.
*   **(x, y)**: Source position (cm).
*   **Param**: Standard deviation (Gaussian) or length/diameter (Flat).

---

## 📊 Output & Visualization

**Console Output:**
*   **Reflectance/Transmittance (Chatterjee/SL)**: R and T values calculated based on different definitions.
*   **MOP**: Mean Optical Path.
*   **User Time**: Simulation duration.

**File Output (.mco):**
Contains grid data for absorption, reflectance, and transmittance.
*   **N_phR / N_phT**: Photon counts on PDs.
*   **summary.csv**: (GPU Version) A concise summary of `Rd` and `Tt` for all runs.

**Visualization (MATLAB):**
1.  Open `get_mop.m`.
2.  Set the filename (e.g., `960.mco`).
3.  Run the script to plot:
    *   Photon weight distribution.
    *   Detector and Light source positions.

---

## 🔧 Programming Notes

*   **Logic**: Core photon transport (Hop/Drop/Spin) is consistent across CPU and GPU versions.
*   **GPU Differences**:
    *   Uses `curand` for parallel random numbers.
    *   Uses atomic operations (`atomicAdd`) for result accumulation.
    *   Flattened arrays for efficient memory access.
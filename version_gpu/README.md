# MCML GPU Version

This directory contains a CUDA-accelerated port of the original MCML (Monte Carlo Multi-Layered) simulation code.

## Key Changes from CPU Version

1.  **CUDA Kernels**: The core photon transport logic (`Hop`, `Drop`, `Spin` etc.) has been moved to CUDA device functions.
2.  **Parallelization**: 
    - Millions of photons are simulated in parallel using thousands of GPU threads.
    - Each thread handles an independent photon.
    - `curand` is used for parallel random number generation (replacing the serial RNG).
3.  **Memory Management**:
    - Input parameters are copied to constant/device memory.
    - Output arrays (e.g., `A_rz`, `Rd_ra`) are **flattened** into 1D arrays to simplify memory allocation and access on the GPU.
    - **Atomic Operations** (`atomicAdd`) are used to safely accumulate results (Absorption, Reflectance, Transmittance) from multiple threads into shared global memory arrays.
4.  **I/O Headers**:
    - The `main` function was updated to correctly parse the header (Version, NumRuns) and loop through all simulation runs defined in the input file.
5.  **Performance**:
    - Significant speedup is achieved by leveraging the massive parallelism of GPUs (e.g., NVIDIA A100).

## Usage

### Prerequisites
- NVIDIA GPU (Compute Capability 5.0+)
- CUDA Toolkit (tested with 11.8)
- GCC (compatible version, e.g., gcc 11)

### Compilation
Use the provided `Makefile`:

```bash
cd version_gpu
make
```

This will generate the `mcml_gpu` executable.

### Running a Simulation

```bash
./mcml_gpu <input_file.mci>
```

Example:
```bash
./mcml_gpu output_5nm.mci
```

### SLURM Submission (for Clusters)
A sample SLURM script `run_gpu.slurm` is provided for running on HPC clusters (like the LIP6 convergence cluster).

```bash
sbatch run_gpu.slurm
```

## Output
The program produces `.mco` output files compatible with the original MCML format, which can be analyzed using standard MATLAB scripts (e.g., `get_mop.m`).

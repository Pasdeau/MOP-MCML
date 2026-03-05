# MOP-MCML: Monte Carlo Multi-Layered Simulation

MOP-MCML extends the classical **MCML** (Multi-Layered Monte Carlo) light transport simulation with two key additions:

1. **Photodiode (PD) Geometry**: Detectors at the reflection and/or transmission side are explicitly defined in the input file. The simulation counts photons actually received by each PD (`NphR`, `NphT`, `NphOP`).
2. **Optical Path (OP) / Mean Optical Path (MOP)**: Spatial distribution of photon paths from source to PD, stored as `OP(z,r)` (2D) or `OP_3D(z,y,x)` (3D), and visualized as a "banana curve".

**Authors:**
- Songlin Li – Original modifications
- Wenzheng Wang – Improvements, GPU port, 3D extension

**Supervisors:**
- Dr. Julien Denoulet (Sorbonne Université, LIP6)
- Prof. Sylvain Feruglio (Sorbonne Université, LIP6)

---

## 📂 Project Structure

```
MOP-MCML/
├── mcmlmain.c / mcmlgo.c / mcmlio.c / mcmlnr.c / mcml.h  ← Version 1: 2D CPU
├── version_gpu/       ← Version 2: 2D GPU (CUDA)
├── version_3d_cpu/    ← Version 3: 3D CPU
├── version_3d_gpu/    ← Version 4: 3D GPU (CUDA)
├── look_mop.m / mop_read_mco.m / mop_plot_mco.m          ← 2D Visualization (MATLAB)
├── look_mop_3d.m / mop_read_mco_3d.m                     ← 3D Visualization (MATLAB)
├── mop_extract_mop.m                                      ← MOP curve extraction
└── compare_mco.py                                         ← CPU/GPU result comparison
```

---

## 🚀 Version 1: 2D CPU (Root Directory)

**Purpose**: Baseline version. Easiest to debug. Use this first to validate your `.mci` parameters, PD/source positions, and confirm the expected "banana-shaped" OP distribution.

**Coordinates**: Cylindrical symmetry (`z`, `r`). Outputs `OP(z,r)`.

### Compilation

**Linux / macOS:**
```bash
gcc -O3 -o mcml mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -lm
./mcml test.mci
```

**Windows (Visual Studio 2022):**
1. Open `projet_MCML.sln`
2. Build Solution (Ctrl+Shift+B)

### PD Mode Flag

| `PD_MODE` | Behavior |
|---|---|
| `1` | Reflectance PD only |
| `2` | Transmittance PD only |
| `3` (default) | Both R + T |

```bash
gcc -O3 -o mcml_Ronly mcmlmain.c mcmlgo.c mcmlio.c mcmlnr.c -lm -DPD_MODE=1
```

---

## ⚡ Version 2: 2D GPU (`version_gpu/`)

**Purpose**: CUDA-accelerated version of the 2D simulation. Same physics and output format as the 2D CPU version. Use for large-scale parameter sweeps.

**Coordinates**: Same as 2D CPU — `OP(z,r)`, output block `OP`.

### Performance

| | CPU Version | GPU Version |
|---|---|---|
| Time / run | ~70 s | **~0.4 s** (A100) |
| Speedup | 1× | **~175×** |
| Transport | Serial | Parallel (CUDA) |

*Benchmark: 2,000,000 photons/run.*

### Compilation

```bash
cd version_gpu
make          # produces ./mcml_gpu
```

**Requirements:** CUDA Toolkit (≥ 11.8), GCC (compatible), NVIDIA GPU (Compute Capability 5.0+).

### Running

```bash
./mcml_gpu test.mci
```

### HPC / SLURM

```bash
sbatch version_gpu/run_gpu.slurm           # single file
sbatch version_gpu/run_batch.slurm         # multiple files sequentially
sbatch version_gpu/run_gpu_variants.slurm  # Ronly / Tonly / RT modes
```

```bash
# Minimal SLURM template
#!/bin/bash
#SBATCH --partition=convergence
#SBATCH --gres=gpu:a100_3g.40gb:1

module purge
module load gcc/11 cuda/11.8
./mcml_gpu your_file.mci
```

---

## 🧊 Version 3: 3D CPU (`version_3d_cpu/`)

**Purpose**: Extends OP from 2D cylindrical to full 3D Cartesian space (`x, y, z`). CPU version is intended for correctness validation and small-scale runs.

**Coordinates**: 3D Cartesian. Outputs `OP_3D(z,y,x)` (block `OP_3D` in `.mco`).

### Input Differences vs. 2D

The `.mci` file adds `dx dy` and `Nx Ny` parameters:

```
dz dr dx dy        # voxel sizes (cm)
Nz Nr Na Nx Ny     # grid dimensions (Na unused for OP, kept for compatibility)
```

See `version_3d_cpu/test_3d.mci` for a working example.

### Compilation & Running

```bash
cd version_3d_cpu
bash run_3d_variants.sh   # compiles and runs Ronly / Tonly / RT modes
```

The script also runs both ASCII (A) and Binary (B) output formats using `test_3d.mci` and `test_3d_B.mci`.

---

## 🚀🧊 Version 4: 3D GPU (`version_3d_gpu/`)

**Purpose**: Full 3D OP simulation with CUDA acceleration. Use for large-scale 3D parameter sweeps where the CPU version is too slow.

**Coordinates**: Same as 3D CPU — `OP_3D(z,y,x)`.

### Running on HPC

```bash
sbatch version_3d_gpu/run_3d_gpu_variants.slurm
# Runs Ronly / Tonly / RT modes, outputs: 3d_gpu_Ronly.mco, etc.
```

**Requirements:** Same as 2D GPU (CUDA ≥ 11.8, NVIDIA GPU).

---

## 📝 Input Format (`.mci` Files)

All 4 versions share the same MCML-style `.mci` format. The last two lines define the PD geometry and light source:

**PD line (2nd to last):**
```
Rx Ry Rl  Tx Ty Tl
```
- `(Rx, Ry)`: Reflectance PD center (cm)
- `Rl`: Reflectance PD side length (cm)
- `(Tx, Ty)`: Transmittance PD center (cm)
- `Tl`: Transmittance PD side length (cm)

**Light source (last line):**
```
Type  x  y  Param
```
- `Type`: `1` = Point, `2` = Gaussian, `3` = Flat
- `(x, y)`: Source position (cm)
- `Param`: σ (Gaussian) or radius/length (Flat), in cm

---

## 📊 Output & Visualization

### `.mco` File Structure

All versions produce `.mco` files in standard MCML format with MOP-MCML additions at the end:

| Block | Versions | Description |
|---|---|---|
| `InParm`, `RAT`, `A_z`, `Rd_r`, `Tt_r`, ... | All | Standard MCML output |
| `OP` | 2D (CPU/GPU) | Optical path distribution `OP(z,r)` |
| `OP_3D` | 3D (CPU/GPU) | 3D optical path `OP_3D(z,y,x)` |
| `NphR`, `NphT`, `NphOP` | All | Photon counts on PDs |

> **Normalization**: The OP arrays are stored as raw sums. Normalize by `NphOP` at the read/plot stage (handled automatically by the MATLAB scripts).

### Visualization (MATLAB)

**2D:**
```matlab
% Edit mco_file at the top, then run:
look_mop        % OP banana plot + optional MOP curve overlay
```

**3D:**
```matlab
% Edit mco_file and SliceFactor at the top, then run:
look_mop_3d     % interactive xyz slice viewer with sliders
```

---

## 🔍 CPU/GPU Comparison Tool

`compare_mco.py` computes Pearson correlation and relative error between two `.mco` files:

```bash
python3 compare_mco.py ref.mco test.mco
# Compares: RAT scalars, Rd_r / Tt_r radial curves, OP / OP_3D arrays
```

---

## 📖 Recommended Workflow

For any new material / PD geometry, follow this order:

1. **2D CPU** — validate `.mci` parameters, PD/source positions, check R/T trends.
2. **2D GPU** — run with identical `.mci` and photon count; compare with `compare_mco.py`.
3. **3D CPU (small run)** — confirm `OP_3D` dimensions and visualization are correct.
4. **3D GPU** — run large-scale batches on HPC; keep fixed module versions for reproducibility.

---

## License

BSD 3-Clause License — see [LICENSE](LICENSE).

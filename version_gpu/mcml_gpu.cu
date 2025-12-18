/***********************************************************
 *  MCML GPU Version (CUDA)
 *  Main logic and Kernels.
 ***********************************************************/

#include "mcml_gpu.h"
#include <cuda_runtime.h>
#include <curand_kernel.h>

#define GPU_CHECK(call)                                                        \
  do {                                                                         \
    cudaError_t err = call;                                                    \
    if (err != cudaSuccess) {                                                  \
      fprintf(stderr, "CUDA Error: %s at %s:%d\n", cudaGetErrorString(err),    \
              __FILE__, __LINE__);                                             \
      exit(1);                                                                 \
    }                                                                          \
  } while (0)

/***********************************************************
 *  Device Constants / Helpers
 ***********************************************************/

__device__ double AtomicAddDouble(double *address, double val) {
  unsigned long long int *address_as_ull = (unsigned long long int *)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
                    __double_as_longlong(val + __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

/***********************************************************
 *  Device Functions corresponding to mcmlgo.c
 ***********************************************************/

__device__ double Rspecular(LayerStruct *Layers) {
  double r1, r2;
  double temp;

  temp = (Layers[0].n - Layers[1].n) / (Layers[0].n + Layers[1].n);
  r1 = temp * temp;

  if ((Layers[1].mua == 0.0) && (Layers[1].mus == 0.0)) { /* glass layer. */
    temp = (Layers[1].n - Layers[2].n) / (Layers[1].n + Layers[2].n);
    r2 = temp * temp;
    r1 = r1 + (1 - r1) * (1 - r1) * r2 / (1 - r1 * r2);
  }
  return r1;
}

__device__ double RFresnel(double n1, double n2, double ca1, double *ca2_Ptr) {
  double r;
  if (n1 == n2) {
    *ca2_Ptr = ca1;
    r = 0.0;
  } else if (ca1 > (1.0 - 1.0E-12)) {
    *ca2_Ptr = ca1;
    r = (n2 - n1) / (n2 + n1);
    r *= r;
  } else if (ca1 < 1.0E-6) {
    *ca2_Ptr = 0.0;
    r = 1.0;
  } else {
    double sa1 = sqrt(1 - ca1 * ca1);
    double sa2 = n1 * sa1 / n2;
    if (sa2 >= 1.0) {
      *ca2_Ptr = 0.0;
      r = 1.0;
    } else {
      double ca2 = sqrt(1 - sa2 * sa2);
      *ca2_Ptr = ca2;
      double cap = ca1 * ca2 - sa1 * sa2;
      double cam = ca1 * ca2 + sa1 * sa2;
      double sap = sa1 * ca2 + ca1 * sa2;
      double sam = sa1 * ca2 - ca1 * sa2;
      r = 0.5 * sam * sam * (cam * cam + cap * cap) / (sap * sap * cam * cam);
    }
  }
  return r;
}

__device__ void SpinTheta(double g, double *cost, curandState *state) {
  if (g == 0.0) {
    *cost = 2.0 * curand_uniform_double(state) - 1.0;
  } else {
    double temp = (1 - g * g) / (1 - g + 2 * g * curand_uniform_double(state));
    *cost = (1 + g * g - temp * temp) / (2 * g);
    if (*cost < -1.0)
      *cost = -1.0;
    else if (*cost > 1.0)
      *cost = 1.0;
  }
}

__device__ void Spin(double g, PhotonStruct *Photon_Ptr, curandState *state) {
  double cost;
  SpinTheta(g, &cost, state);
  double sint = sqrt(1.0 - cost * cost);

  double psi = 2.0 * PI * curand_uniform_double(state);
  double cosp = cos(psi);
  double sinp = (psi < PI) ? sqrt(1.0 - cosp * cosp) : -sqrt(1.0 - cosp * cosp);

  double ux = Photon_Ptr->ux;
  double uy = Photon_Ptr->uy;
  double uz = Photon_Ptr->uz;

  if (fabs(uz) > (1.0 - 1.0E-12)) {
    Photon_Ptr->ux = sint * cosp;
    Photon_Ptr->uy = sint * sinp;
    Photon_Ptr->uz = cost * SIGN(uz);
  } else {
    double temp = sqrt(1.0 - uz * uz);
    Photon_Ptr->ux = sint * (ux * uz * cosp - uy * sinp) / temp + ux * cost;
    Photon_Ptr->uy = sint * (uy * uz * cosp + ux * sinp) / temp + uy * cost;
    Photon_Ptr->uz = -sint * cosp * temp + uz * cost;
  }
}

__device__ void Hop(PhotonStruct *Photon_Ptr) {
  double s = Photon_Ptr->s;
  Photon_Ptr->x += s * Photon_Ptr->ux;
  Photon_Ptr->y += s * Photon_Ptr->uy;
  Photon_Ptr->z += s * Photon_Ptr->uz;
  Photon_Ptr->OP += s;
}

__device__ void StepSizeInGlass(PhotonStruct *Photon_Ptr, InputStruct *In_Ptr,
                                LayerStruct *Layers) {
  double dl_b;
  short layer = Photon_Ptr->layer;
  double uz = Photon_Ptr->uz;

  if (uz > 0.0)
    dl_b = (Layers[layer].z1 - Photon_Ptr->z) / uz;
  else if (uz < 0.0)
    dl_b = (Layers[layer].z0 - Photon_Ptr->z) / uz;
  else
    dl_b = 0.0;

  Photon_Ptr->s = dl_b;
}

__device__ void StepSizeInTissue(PhotonStruct *Photon_Ptr, InputStruct *In_Ptr,
                                 LayerStruct *Layers, curandState *state) {
  short layer = Photon_Ptr->layer;
  double mua = Layers[layer].mua;
  double mus = Layers[layer].mus;

  if (Photon_Ptr->sleft == 0.0) {
    double rnd;
    do
      rnd = curand_uniform_double(state);
    while (rnd <= 0.0);
    Photon_Ptr->s = -log(rnd) / (mua + mus);
  } else {
    Photon_Ptr->s = Photon_Ptr->sleft / (mua + mus);
    Photon_Ptr->sleft = 0.0;
  }
}

__device__ bool HitBoundary(PhotonStruct *Photon_Ptr, InputStruct *In_Ptr,
                            LayerStruct *Layers) {
  double dl_b = 0.0;
  short layer = Photon_Ptr->layer;
  double uz = Photon_Ptr->uz;

  if (uz > 0.0)
    dl_b = (Layers[layer].z1 - Photon_Ptr->z) / uz;
  else if (uz < 0.0)
    dl_b = (Layers[layer].z0 - Photon_Ptr->z) / uz;

  if (uz != 0.0 && Photon_Ptr->s > dl_b) {
    double mut = Layers[layer].mua + Layers[layer].mus;
    Photon_Ptr->sleft = (Photon_Ptr->s - dl_b) * mut;
    Photon_Ptr->s = dl_b;
    return true;
  }
  return false;
}

__device__ void Drop(InputStruct *In_Ptr, LayerStruct *Layers,
                     PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr) {
  double x = Photon_Ptr->x;
  double y = Photon_Ptr->y;
  short iz = (short)(Photon_Ptr->z / In_Ptr->dz);
  if (iz > In_Ptr->nz - 1)
    iz = In_Ptr->nz - 1;

  short ir = (short)(sqrt(x * x + y * y) / In_Ptr->dr);
  if (ir > In_Ptr->nr - 1 || ir < 0)
    ir = In_Ptr->nr - 1;

  short layer = Photon_Ptr->layer;
  double mua = Layers[layer].mua;
  double mus = Layers[layer].mus;
  double dwa = Photon_Ptr->w * mua / (mua + mus);
  Photon_Ptr->w -= dwa;

  AtomicAddDouble(&Out_Ptr->A_rz[ir * In_Ptr->nz + iz], dwa);

  // We don't record step path (cpt, sOut, etc) here as that was for
  // debug/special purposes and requires huge memory
}

__device__ void Roulette(PhotonStruct *Photon_Ptr, InputStruct *In_Ptr,
                         curandState *state) {
  if (Photon_Ptr->w == 0.0) {
    Photon_Ptr->dead = true;
  } else if (curand_uniform_double(state) < CHANCE) {
    Photon_Ptr->w /= CHANCE;
  } else {
    Photon_Ptr->dead = true;
  }
}

__device__ void RecordR(double Refl, InputStruct *In_Ptr,
                        PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr) {
  double x = Photon_Ptr->x, y = Photon_Ptr->y;
  short ir = (short)(sqrt(x * x + y * y) / In_Ptr->dr);
  if (ir > In_Ptr->nr - 1 || ir < 0)
    ir = In_Ptr->nr - 1;

  short ia = (short)(acos(-Photon_Ptr->uz) / In_Ptr->da);
  if (ia > In_Ptr->na - 1)
    ia = In_Ptr->na - 1;

  AtomicAddDouble(&Out_Ptr->Rd_ra[ir * In_Ptr->na + ia],
                  Photon_Ptr->w * (1.0 - Refl));
  Photon_Ptr->w *= Refl;
}

__device__ void RecordT(double Refl, InputStruct *In_Ptr,
                        PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr) {
  double x = Photon_Ptr->x, y = Photon_Ptr->y;
  short ir = (short)(sqrt(x * x + y * y) / In_Ptr->dr);
  if (ir > In_Ptr->nr - 1 || ir < 0)
    ir = In_Ptr->nr - 1;

  short ia = (short)(acos(Photon_Ptr->uz) / In_Ptr->da);
  if (ia > In_Ptr->na - 1)
    ia = In_Ptr->na - 1;

  AtomicAddDouble(&Out_Ptr->Tt_ra[ir * In_Ptr->na + ia],
                  Photon_Ptr->w * (1.0 - Refl));
  Photon_Ptr->w *= Refl;
}

__device__ void CrossUpOrNot(InputStruct *In_Ptr, LayerStruct *Layers,
                             PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr,
                             curandState *state) {
  double uz = Photon_Ptr->uz;
  double uz1;
  double r = 0.0;
  short layer = Photon_Ptr->layer;
  double ni = Layers[layer].n;
  double nt = Layers[layer - 1].n;

  if (-uz <= Layers[layer].cos_crit0)
    r = 1.0;
  else
    r = RFresnel(ni, nt, -uz, &uz1);

  if (curand_uniform_double(state) > r) {
    if (layer == 1) {
      Photon_Ptr->uz = -uz1;
      RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
      Photon_Ptr->dead = true;
    } else {
      Photon_Ptr->layer--;
      Photon_Ptr->ux *= ni / nt;
      Photon_Ptr->uy *= ni / nt;
      Photon_Ptr->uz = -uz1;
    }
  } else {
    Photon_Ptr->uz = -uz;
  }
}

__device__ void CrossDnOrNot(InputStruct *In_Ptr, LayerStruct *Layers,
                             PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr,
                             curandState *state) {
  double uz = Photon_Ptr->uz;
  double uz1;
  double r = 0.0;
  short layer = Photon_Ptr->layer;
  double ni = Layers[layer].n;
  double nt = Layers[layer + 1].n;

  if (uz <= Layers[layer].cos_crit1)
    r = 1.0;
  else
    r = RFresnel(ni, nt, uz, &uz1);

  if (curand_uniform_double(state) > r) {
    if (layer == In_Ptr->num_layers) {
      Photon_Ptr->uz = uz1;
      RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
      Photon_Ptr->dead = true;
    } else {
      Photon_Ptr->layer++;
      Photon_Ptr->ux *= ni / nt;
      Photon_Ptr->uy *= ni / nt;
      Photon_Ptr->uz = uz1;
    }
  } else {
    Photon_Ptr->uz = -uz;
  }
}

__device__ void CrossOrNot(InputStruct *In_Ptr, LayerStruct *Layers,
                           PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr,
                           curandState *state) {
  if (Photon_Ptr->uz < 0.0)
    CrossUpOrNot(In_Ptr, Layers, Photon_Ptr, Out_Ptr, state);
  else
    CrossDnOrNot(In_Ptr, Layers, Photon_Ptr, Out_Ptr, state);
}

__device__ void HopInGlass(InputStruct *In_Ptr, LayerStruct *Layers,
                           PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr,
                           curandState *state) {
  if (Photon_Ptr->uz == 0.0) {
    Photon_Ptr->dead = true;
  } else {
    StepSizeInGlass(Photon_Ptr, In_Ptr, Layers);
    Hop(Photon_Ptr);
    CrossOrNot(In_Ptr, Layers, Photon_Ptr, Out_Ptr, state);
  }
}

__device__ void HopDropSpinInTissue(InputStruct *In_Ptr, LayerStruct *Layers,
                                    PhotonStruct *Photon_Ptr,
                                    OutStruct *Out_Ptr, curandState *state) {
  StepSizeInTissue(Photon_Ptr, In_Ptr, Layers, state);

  if (HitBoundary(Photon_Ptr, In_Ptr, Layers)) {
    Hop(Photon_Ptr);
    CrossOrNot(In_Ptr, Layers, Photon_Ptr, Out_Ptr, state);
  } else {
    Hop(Photon_Ptr);
    Drop(In_Ptr, Layers, Photon_Ptr, Out_Ptr);
    Spin(Layers[Photon_Ptr->layer].g, Photon_Ptr, state);
  }
}

__device__ void HopDropSpin(InputStruct *In_Ptr, LayerStruct *Layers,
                            PhotonStruct *Photon_Ptr, OutStruct *Out_Ptr,
                            curandState *state) {
  short layer = Photon_Ptr->layer;

  if ((Layers[layer].mua == 0.0) && (Layers[layer].mus == 0.0))
    HopInGlass(In_Ptr, Layers, Photon_Ptr, Out_Ptr, state);
  else
    HopDropSpinInTissue(In_Ptr, Layers, Photon_Ptr, Out_Ptr, state);

  if (Photon_Ptr->w < In_Ptr->Wth && !Photon_Ptr->dead)
    Roulette(Photon_Ptr, In_Ptr, state);
}

__device__ void LaunchPhoton(LayerStruct *Layers, PhotonStruct *Photon_Ptr,
                             short lightType, double light_x, double light_y,
                             double light_l, curandState *state) {
  Photon_Ptr->w = 1.0 - Rspecular(Layers);
  Photon_Ptr->dead = false;
  Photon_Ptr->layer = 1;
  Photon_Ptr->s = 0.0;
  Photon_Ptr->sleft = 0.0;

  double A, B, C, rx, ry;
  double mu = 0.0, sigma = light_l;

  if (lightType == 2) {
    // Box-Muller
    A = sqrt((-2.0) * log(curand_uniform_double(state)));
    B = 2.0 * PI * curand_uniform_double(state);
    C = A * cos(B);
    rx = mu + C * sigma;

    A = sqrt((-2.0) * log(curand_uniform_double(state)));
    B = 2.0 * PI * curand_uniform_double(state);
    C = A * cos(B);
    ry = mu + C * sigma;

    Photon_Ptr->x = light_x + rx;
    Photon_Ptr->y = light_y + ry;
  } else if (lightType == 3) {
    double half = 0.5 * sigma;
    Photon_Ptr->x =
        light_x + (curand_uniform_double(state) * 2.0 * half) - half;
    Photon_Ptr->y =
        light_y + (curand_uniform_double(state) * 2.0 * half) - half;
  } else {
    Photon_Ptr->x = light_x;
    Photon_Ptr->y = light_y;
  }

  Photon_Ptr->z = 0.0;
  Photon_Ptr->ux = 0.0;
  Photon_Ptr->uy = 0.0;
  Photon_Ptr->uz = 1.0;

  if ((Layers[1].mua == 0.0) && (Layers[1].mus == 0.0)) {
    Photon_Ptr->layer = 2;
    Photon_Ptr->z = Layers[2].z0;
  }
}

/***********************************************************
 *  RNG Setup Kernel
 ***********************************************************/
__global__ void SetupRNG(curandState *states, unsigned long long seed,
                         long num_photons) {
  long idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx < num_photons) {
    curand_init(seed, idx, 0, &states[idx]);
  }
}

/***********************************************************
 *  Main MCML Kernel
 ***********************************************************/
__global__ void MCML_Kernel(InputStruct In_Ptr, LayerStruct *Layers,
                            OutStruct Out_Ptr, curandState *states) {
  long idx = blockIdx.x * blockDim.x + threadIdx.x;
  if (idx >= In_Ptr.num_photons)
    return;

  curandState localState = states[idx];
  PhotonStruct photon;

  LaunchPhoton(Layers, &photon, In_Ptr.lightType, In_Ptr.light_x,
               In_Ptr.light_y, In_Ptr.light_l, &localState);

  while (!photon.dead) {
    HopDropSpin(&In_Ptr, Layers, &photon, &Out_Ptr, &localState);
  }

  states[idx] = localState; // Save state if we were to continue (not reused
                            // here but good practice)
}

/***********************************************************
 *  Main
 ***********************************************************/
int main(int argc, char *argv[]) {
  char inputs_fname[STRLEN] = "";
  FILE *file;
  InputStruct In_Parm;
  OutStruct Out_Parm;

  if (argc >= 2)
    strcpy(inputs_fname, argv[1]);
  file = GetFile(inputs_fname);

  CheckParm(file, &In_Parm);          // Logic check only
  short num_runs = ReadNumRuns(file); // Skip header (Version and NumRuns)

  // Open summary CSV (overwrite)
  const char *summary_name = "summary.csv";
  FILE *summary_fp = fopen(summary_name, "w");
  if (!summary_fp) {
    fprintf(
        stderr,
        "Warning: cannot open %s for writing. Summary CSV will be skipped.\n",
        summary_name);
  } else {
    // Header
    fprintf(summary_fp, "output,Rd,Tt\n");
    fflush(summary_fp);
  }

  for (short i_run = 1; i_run <= num_runs; i_run++) {
    printf("\nProcessing Run %hd of %hd\n", i_run, num_runs);
    ReadParm(file, &In_Parm); // Actually read run data
    InitOutputData(In_Parm, &Out_Parm);

    printf("Simulating %ld photons on GPU...\n", In_Parm.num_photons);

    // Device Memory Allocation
    LayerStruct *d_Layers;
    GPU_CHECK(
        cudaMalloc(&d_Layers, (In_Parm.num_layers + 2) * sizeof(LayerStruct)));
    GPU_CHECK(cudaMemcpy(d_Layers, In_Parm.layerspecs,
                         (In_Parm.num_layers + 2) * sizeof(LayerStruct),
                         cudaMemcpyHostToDevice));

    // Output Memory Allocation
    OutStruct d_Out_Parm;
    // Flattened arrays
    GPU_CHECK(cudaMalloc(&(d_Out_Parm.Rd_ra),
                         In_Parm.nr * In_Parm.na * sizeof(double)));
    GPU_CHECK(cudaMemset(d_Out_Parm.Rd_ra, 0,
                         In_Parm.nr * In_Parm.na * sizeof(double)));

    GPU_CHECK(cudaMalloc(&(d_Out_Parm.A_rz),
                         In_Parm.nr * In_Parm.nz * sizeof(double)));
    GPU_CHECK(cudaMemset(d_Out_Parm.A_rz, 0,
                         In_Parm.nr * In_Parm.nz * sizeof(double)));

    GPU_CHECK(cudaMalloc(&(d_Out_Parm.Tt_ra),
                         In_Parm.nr * In_Parm.na * sizeof(double)));
    GPU_CHECK(cudaMemset(d_Out_Parm.Tt_ra, 0,
                         In_Parm.nr * In_Parm.na * sizeof(double)));

    // RNG States
    curandState *d_states;
    GPU_CHECK(cudaMalloc(&d_states, In_Parm.num_photons * sizeof(curandState)));

    // Launch Config
    int blockSize = 256;
    int numBlocks = (In_Parm.num_photons + blockSize - 1) / blockSize;

    // init RNG
    SetupRNG<<<numBlocks, blockSize>>>(
        d_states, time(NULL) + i_run,
        In_Parm.num_photons); // Seed varies by run
    GPU_CHECK(cudaDeviceSynchronize());

    clock_t start = clock();

    // Run Simulation
    MCML_Kernel<<<numBlocks, blockSize>>>(In_Parm, d_Layers, d_Out_Parm,
                                          d_states);
    GPU_CHECK(cudaDeviceSynchronize());

    clock_t end = clock();
    double duration = (double)(end - start) / CLOCKS_PER_SEC;
    printf("Simulation done in %.2f sec\n", duration);

    // Copy back results
    GPU_CHECK(cudaMemcpy(Out_Parm.Rd_ra, d_Out_Parm.Rd_ra,
                         In_Parm.nr * In_Parm.na * sizeof(double),
                         cudaMemcpyDeviceToHost));
    GPU_CHECK(cudaMemcpy(Out_Parm.A_rz, d_Out_Parm.A_rz,
                         In_Parm.nr * In_Parm.nz * sizeof(double),
                         cudaMemcpyDeviceToHost));
    GPU_CHECK(cudaMemcpy(Out_Parm.Tt_ra, d_Out_Parm.Tt_ra,
                         In_Parm.nr * In_Parm.na * sizeof(double),
                         cudaMemcpyDeviceToHost));

    // Post-processing on Host
    SumScaleResult(In_Parm, &Out_Parm);
    WriteResult(In_Parm, Out_Parm, (char *)"GPU Simulation");

    // Append to summary.csv and print a line to terminal
    if (summary_fp) {
      // CSV: output,Rd,Tt
      fprintf(summary_fp, "%s,%.10g,%.10g\n", In_Parm.out_fname, Out_Parm.Rd,
              Out_Parm.Tt);
      fflush(summary_fp);
    }
    printf("[SUMMARY] %s -> Rd=%.6g, Tt=%.6g\n", In_Parm.out_fname, Out_Parm.Rd,
           Out_Parm.Tt);

    // Cleanup
    cudaFree(d_Layers);
    cudaFree(d_Out_Parm.Rd_ra);
    cudaFree(d_Out_Parm.A_rz);
    cudaFree(d_Out_Parm.Tt_ra);
    cudaFree(d_states);
    FreeData(In_Parm, &Out_Parm);
  }

  if (summary_fp) {
    fclose(summary_fp);
    printf("All runs done. Summary written to %s\n", summary_name);
  }

  fclose(file);
  return 0;
}

#ifndef MCML_GPU_H
#define MCML_GPU_H

#include <ctype.h>
#include <math.h>
#include <stdbool.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define PI 3.1415926
#define WEIGHT 1E-4
#define CHANCE 0.1
#define STRLEN 256
#define SIGN(x) ((x) >= 0 ? 1 : -1)

typedef bool Boolean;

/****************** Structures *****************************/

typedef struct {
  double x, y, z;
  double ux, uy, uz;
  double w;
  bool dead;
  short layer;
  double s;
  double sleft;
  double OP;
} PhotonStruct;

typedef struct {
  double z0, z1;
  double n;
  double mua;
  double mus;
  double g;
  double cos_crit0, cos_crit1;
} LayerStruct;

typedef struct {
  char out_fname[STRLEN];
  char out_fformat;
  long num_photons;
  double Wth;

  double dz;
  double dr;
  double da;
  short nz;
  short nr;
  short na;

  short num_layers;
  LayerStruct *layerspecs; /* Host pointer, will be copied to constant memory or
                              device memory */

  double PD_Rx;
  double PD_Ry;
  double PD_Tx;
  double PD_Ty;
  double PD_Tz;
  double PD_Rl;
  double PD_Tl;

  short lightType;
  double light_x;
  double light_y;
  double light_l;
} InputStruct;

/* Flattened Output Structure for GPU */
typedef struct {
  double Rsp;

  /* Flattened arrays: size = nr * na */
  double *Rd_ra;

  double *Rd_r; /* size = nr */
  double *Rd_a; /* size = na */
  double Rd;

  /* Flattened arrays: size = nr * nz */
  double *A_rz;

  double *A_z; /* size = nz */
  double *A_l; /* size = num_layers + 2 */
  double A;

  /* Flattened arrays: size = nr * na */
  double *Tt_ra;

  double *Tt_r; /* size = nr */
  double *Tt_a; /* size = na */
  double Tt;

  /* Flattened arrays: size = nr * nz */
  double *OP;

  int photonsnbrR;
  int photonsnbrT;
  int nbrphotons;
} OutStruct;

/* Function Prototypes for IO */
#ifdef __cplusplus
extern "C" {
#endif

void ShowVersion(const char *version);
FILE *GetFile(char *Fname);
short ReadNumRuns(FILE *File_Ptr);
void CheckParm(FILE *File_Ptr, InputStruct *In_Ptr);
void ReadParm(FILE *File_Ptr, InputStruct *In_Ptr);
void InitOutputData(InputStruct In_Parm, OutStruct *Out_Ptr);
void FreeData(InputStruct In_Parm, OutStruct *Out_Ptr);
void SumScaleResult(InputStruct In_Parm, OutStruct *Out_Ptr);
void WriteResult(InputStruct In_Parm, OutStruct Out_Parm, char *TimeReport);
void nrerror(char *error_text);

#ifdef __cplusplus
}
#endif

#endif

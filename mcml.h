#pragma once
/***********************************************************
 *  Monte Carlo simulation of photon distribution in
 *  multi-layered turbid media in ANSI Standard C.
 *
 *  Minimal refinements for macOS/Clang (C11):
 *  - Use <stdbool.h> and typedef bool Boolean
 *  - Provide prototype for ShowVersion
 *  - Keep original structures and constants
 ****/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <time.h>
#include <string.h>
#include <ctype.h>
#include <stdbool.h>

#define PI 3.1415926
#define WEIGHT 1E-4     /* Critical weight for roulette. */
#define CHANCE 0.1      /* Chance of roulette survival. */
#define STRLEN 256      /* String length. */

typedef bool Boolean;

#define SIGN(x) ((x)>=0 ? 1:-1)

/****************** Structures *****************************/

/****
 * Structure used to describe a photon packet.
 ****/
typedef struct {
    double x, y, z;     /* Cartesian coordinates.[cm] */
    double ux, uy, uz;  /* directional cosines of a photon. */
    double w;           /* weight. */
    Boolean dead;       /* 1 if photon is terminated. */
    short layer;        /* index to layer where the photon packet resides. */
    double s;           /* current step size. [cm]. */
    double sleft;       /* step size left. dimensionless [-]. */
    double OP;          /* optical path accumulator. */
} PhotonStruct;

/****
 * Structure for a layer.
 ****/
typedef struct {
    double z0, z1;  /* z coordinates of a layer. [cm] */
    double n;       /* refractive index of a layer. */
    double mua;     /* absorption coefficient. [1/cm] */
    double mus;     /* scattering coefficient. [1/cm] */
    double g;       /* anisotropy. */

    double cos_crit0, cos_crit1;
} LayerStruct;

/****
 * Input parameters for each run.
 ****/
typedef struct {
    char   out_fname[STRLEN]; /* output file name. */
    char   out_fformat;       /* 'A' for ASCII, 'B' for binary. */
    long   num_photons;       /* photons to trace. */
    double Wth;               /* roulette threshold. */

    double dz;                /* z grid separation.[cm] */
    double dr;                /* r grid separation.[cm] */
    double da;                /* alpha grid separation. [radian] */
    short  nz;                /* array range 0..nz-1. */
    short  nr;                /* array range 0..nr-1. */
    short  na;                /* array range 0..na-1. */

    short  num_layers;        /* number of layers. */
    LayerStruct* layerspecs;  /* layer parameters. */

    /* Photodetector areas (square windows) */
    double PD_Rx;
    double PD_Ry;
    double PD_Tx;
    double PD_Ty;
    double PD_Tz;
    double PD_Rl;
    double PD_Tl;

    /* Light source parameters */
    short  lightType;   /* 1=point, 2=Gaussian(sigma=light_l), 3=flat(top-hat side=light_l) */
    double light_x;
    double light_y;
    double light_l;
} InputStruct;

/****
 * Output scoring.
 ****/
typedef struct {
    double    Rsp;     /* specular reflectance. [-] */
    double** Rd_ra;   /* 2D diffuse reflectance. [1/(cm2 sr)] */
    double* Rd_r;    /* 1D radial diffuse reflectance. [1/cm2] */
    double* Rd_a;    /* 1D angular diffuse reflectance. [1/sr] */
    double    Rd;      /* total diffuse reflectance. [-] */

    double** A_rz;    /* 2D absorbed weight density over (r,z). [1/cm3] */
    double* A_z;     /* 1D absorption over z. [1/cm] */
    double* A_l;     /* per-layer absorption probability. [-] */
    double    A;       /* total absorption probability. [-] */

    double** Tt_ra;   /* 2D total transmittance. [1/(cm2 sr)] */
    double* Tt_r;    /* 1D radial transmittance. [1/cm2] */
    double* Tt_a;    /* 1D angular transmittance. [1/sr] */
    double    Tt;      /* total transmittance. [-] */

    double** OP;      /* optical path map (example accumulation). */

    int photonsnbrR;
    int photonsnbrT;
    int nbrphotons;
} OutStruct;

/***********************************************************
 * Dynamic memory prototypes (Numerical Recipes style).
 ****/
double* AllocVector(short, short);
double** AllocMatrix(short, short, short, short);
void     FreeVector(double*, short, short);
void     FreeMatrix(double**, short, short, short, short);
void     nrerror(char*);

/* Provide ShowVersion prototype to avoid implicit decl under C11 */
void ShowVersion(const char* version);

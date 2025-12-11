/***********************************************************
 *  Launch, move, and record photon weight.
 *  Refinements:
 *   - unify RNG to ran3/RandomNum (no rand/srand, no globals)
 *   - remove UTF-8 accented comments
 *   - avoid shadow/unused warnings
 ****/

#include "mcml.h"
#include <time.h>
#include <math.h>

 /***********************************************************
  *  A random number generator from Numerical Recipes in C.
  ****/
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC 1.0E-9

static float ran3(int* idum)
{
    static int inext, inextp;
    static long ma[56];
    static int iff = 0;
    long mj, mk;
    int i, ii, k;

    if (*idum < 0 || iff == 0) {
        iff = 1;
        mj = MSEED - (*idum < 0 ? -*idum : *idum);
        mj %= MBIG;
        ma[55] = mj;
        mk = 1;
        for (i = 1; i <= 54; i++) {
            ii = (21 * i) % 55;
            ma[ii] = mk;
            mk = mj - mk;
            if (mk < MZ) mk += MBIG;
            mj = ma[ii];
        }
        for (k = 1; k <= 4; k++)
            for (i = 1; i <= 55; i++) {
                ma[i] -= ma[1 + (i + 30) % 55];
                if (ma[i] < MZ) ma[i] += MBIG;
            }
        inext = 0;
        inextp = 31;
        *idum = 1;
    }
    if (++inext == 56) inext = 1;
    if (++inextp == 56) inextp = 1;
    mj = ma[inext] - ma[inextp];
    if (mj < MZ) mj += MBIG;
    ma[inext] = mj;
    return mj * FAC;
}

#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/***********************************************************
 *  Generate a random number in (0,1).
 ****/
double RandomNum(void)
{
    static Boolean first_time = true;
    static int idum;

    if (first_time) {
#if 0 /* Deterministic test seed (optional) */
        idum = -1;
#else
        idum = -(int)time(NULL) % (1 << 15);
#endif
        (void)ran3(&idum);
        first_time = false;
        idum = 1;
    }
    return (double)ran3(&idum);
}

/* Helpers to unify uniform draws */
static inline double uniform01(void) { /* (0,1) */
    double u;
    do { u = RandomNum(); } while (u <= 0.0);
    return u;
}

/* Generate two independent U(0,1) numbers */
static void UNIFORM(double* p) {
    p[0] = uniform01();
    p[1] = uniform01();
}

/* Uniform over [-length/2, +length/2) */
static double UNIFORM2(double length) {
    double half = 0.5 * length;
    return (uniform01() * (2.0 * half)) - half;
}

/***********************************************************
 *  Compute the specular reflection.
 ****/
double Rspecular(LayerStruct* Layerspecs_Ptr)
{
    double r1, r2;
    double temp;

    temp = (Layerspecs_Ptr[0].n - Layerspecs_Ptr[1].n)
        / (Layerspecs_Ptr[0].n + Layerspecs_Ptr[1].n);
    r1 = temp * temp;

    if ((Layerspecs_Ptr[1].mua == 0.0)
        && (Layerspecs_Ptr[1].mus == 0.0)) { /* glass layer. */
        temp = (Layerspecs_Ptr[1].n - Layerspecs_Ptr[2].n)
            / (Layerspecs_Ptr[1].n + Layerspecs_Ptr[2].n);
        r2 = temp * temp;
        r1 = r1 + (1 - r1) * (1 - r1) * r2 / (1 - r1 * r2);
    }
    return r1;
}

/***********************************************************
 *  Initialize a photon packet and sample launch position.
 ****/
void LaunchPhoton(double Rspecular,
    LayerStruct* Layerspecs_Ptr,
    PhotonStruct* Photon_Ptr,
    short lightType,
    double light_x,
    double light_y,
    double light_l)
{
    Photon_Ptr->w = 1.0 - Rspecular;
    Photon_Ptr->dead = false;
    Photon_Ptr->layer = 1;
    Photon_Ptr->s = 0.0;
    Photon_Ptr->sleft = 0.0;

    double A, B, C, rx, ry;
    double mu = 0.0, sigma = light_l;  /* Gaussian: mean & stddev */
    double uni[2];

    /* Box-Muller for Gaussian coordinates */
    UNIFORM(uni);
    A = sqrt((-2.0) * log(uni[0]));
    B = 2.0 * PI * uni[1];
    C = A * cos(B);
    rx = mu + C * sigma;

    UNIFORM(uni);
    A = sqrt((-2.0) * log(uni[0]));
    B = 2.0 * PI * uni[1];
    C = A * cos(B);
    ry = mu + C * sigma;

    if (lightType == 2) {          /* Gaussian spot */
        Photon_Ptr->x = light_x + rx;
        Photon_Ptr->y = light_y + ry;
    }
    else if (lightType == 3) {   /* flat top-hat (square side=sigma) */
        Photon_Ptr->x = light_x + UNIFORM2(sigma);
        Photon_Ptr->y = light_y + UNIFORM2(sigma);
    }
    else {                       /* point source */
        Photon_Ptr->x = light_x;
        Photon_Ptr->y = light_y;
    }

    Photon_Ptr->z = 0.0;
    Photon_Ptr->ux = 0.0;
    Photon_Ptr->uy = 0.0;
    Photon_Ptr->uz = 1.0;

    if ((Layerspecs_Ptr[1].mua == 0.0) && (Layerspecs_Ptr[1].mus == 0.0)) {
        Photon_Ptr->layer = 2;
        Photon_Ptr->z = Layerspecs_Ptr[2].z0;
    }
}

/***********************************************************
 *  Sample polar deflection angle theta given anisotropy g.
 ****/
double SpinTheta(double g)
{
    double cost;
    if (g == 0.0) {
        cost = 2.0 * RandomNum() - 1.0;
    }
    else {
        double temp = (1 - g * g) / (1 - g + 2 * g * RandomNum());
        cost = (1 + g * g - temp * temp) / (2 * g);
        if (cost < -1.0) cost = -1.0;
        else if (cost > 1.0) cost = 1.0;
    }
    return cost;
}

/***********************************************************
 *  Spin azimuthal angle psi and update direction.
 ****/
void Spin(double g, PhotonStruct* Photon_Ptr)
{
    double cost = SpinTheta(g);
    double sint = sqrt(1.0 - cost * cost);

    double psi = 2.0 * PI * RandomNum();
    double cosp = cos(psi);
    double sinp = (psi < PI) ? sqrt(1.0 - cosp * cosp)
        : -sqrt(1.0 - cosp * cosp);

    double ux = Photon_Ptr->ux;
    double uy = Photon_Ptr->uy;
    double uz = Photon_Ptr->uz;

    if (fabs(uz) > (1.0 - 1.0E-12)) { /* near normal */
        Photon_Ptr->ux = sint * cosp;
        Photon_Ptr->uy = sint * sinp;
        Photon_Ptr->uz = cost * SIGN(uz);
    }
    else {
        double temp = sqrt(1.0 - uz * uz);
        Photon_Ptr->ux = sint * (ux * uz * cosp - uy * sinp) / temp + ux * cost;
        Photon_Ptr->uy = sint * (uy * uz * cosp + ux * sinp) / temp + uy * cost;
        Photon_Ptr->uz = -sint * cosp * temp + uz * cost;
    }
}

/***********************************************************
 *  Move the photon s away in current layer.
 ****/
void Hop(PhotonStruct* Photon_Ptr)
{
    double s = Photon_Ptr->s;
    Photon_Ptr->x += s * Photon_Ptr->ux;
    Photon_Ptr->y += s * Photon_Ptr->uy;
    Photon_Ptr->z += s * Photon_Ptr->uz;
    Photon_Ptr->OP += s; /* accumulate optical path */
}

/***********************************************************
 *  Return step size to boundary in glass (if uz != 0).
 ****/
void StepSizeInGlass(PhotonStruct* Photon_Ptr, InputStruct* In_Ptr)
{
    double dl_b;
    short  layer = Photon_Ptr->layer;
    double uz = Photon_Ptr->uz;

    if (uz > 0.0)
        dl_b = (In_Ptr->layerspecs[layer].z1 - Photon_Ptr->z) / uz;
    else if (uz < 0.0)
        dl_b = (In_Ptr->layerspecs[layer].z0 - Photon_Ptr->z) / uz;
    else
        dl_b = 0.0;

    Photon_Ptr->s = dl_b;
}

/***********************************************************
 *  Pick a step size in tissue.
 ****/
void StepSizeInTissue(PhotonStruct* Photon_Ptr,
    InputStruct* In_Ptr,
    OutStruct* Out_Ptr)
{
    (void)Out_Ptr; /* currently unused here */

    short  layer = Photon_Ptr->layer;
    double mua = In_Ptr->layerspecs[layer].mua;
    double mus = In_Ptr->layerspecs[layer].mus;

    if (Photon_Ptr->sleft == 0.0) {
        double rnd;
        do rnd = RandomNum(); while (rnd <= 0.0);
        Photon_Ptr->s = -log(rnd) / (mua + mus);
    }
    else {
        Photon_Ptr->s = Photon_Ptr->sleft / (mua + mus);
        Photon_Ptr->sleft = 0.0;
    }
}

/***********************************************************
 *  Check if the step will hit a boundary.
 ****/
Boolean HitBoundary(PhotonStruct* Photon_Ptr,
    InputStruct* In_Ptr, OutStruct* Out_Ptr)
{
    (void)Out_Ptr; /* not used here */

    double dl_b = 0.0;
    short  layer = Photon_Ptr->layer;
    double uz = Photon_Ptr->uz;

    if (uz > 0.0)
        dl_b = (In_Ptr->layerspecs[layer].z1 - Photon_Ptr->z) / uz;
    else if (uz < 0.0)
        dl_b = (In_Ptr->layerspecs[layer].z0 - Photon_Ptr->z) / uz;

    if (uz != 0.0 && Photon_Ptr->s > dl_b) {
        double mut = In_Ptr->layerspecs[layer].mua + In_Ptr->layerspecs[layer].mus;
        Photon_Ptr->sleft = (Photon_Ptr->s - dl_b) * mut;
        Photon_Ptr->s = dl_b;
        return true;
    }
    return false;
}

/***********************************************************
 *  Drop photon weight inside tissue.
 ****/
void Drop(InputStruct* In_Ptr,
    PhotonStruct* Photon_Ptr,
    OutStruct* Out_Ptr,
    int* cpt,
    double* sOut,
    short* irOut,
    short* izOut)
{
    (void)In_Ptr;

    double x = Photon_Ptr->x;
    double y = Photon_Ptr->y;
    short  iz = (short)(Photon_Ptr->z / In_Ptr->dz);
    if (iz > In_Ptr->nz - 1) iz = In_Ptr->nz - 1;

    short  ir = (short)(sqrt(x * x + y * y) / In_Ptr->dr);
    if (ir > In_Ptr->nr - 1 || ir < 0) ir = In_Ptr->nr - 1;

    short  layer = Photon_Ptr->layer;
    double mua = In_Ptr->layerspecs[layer].mua;
    double mus = In_Ptr->layerspecs[layer].mus;
    double dwa = Photon_Ptr->w * mua / (mua + mus);
    Photon_Ptr->w -= dwa;

    Out_Ptr->A_rz[ir][iz] += dwa;

    /* Record some step info into tracing arrays (as original did) */
    int idx = *cpt;
    irOut[idx] = ir;
    izOut[idx] = iz;
    sOut[idx] = Photon_Ptr->w;
    *cpt = idx + 1;
}

/***********************************************************
 *  Roulette for low-weight photons.
 ****/
void Roulette(PhotonStruct* Photon_Ptr)
{
    if (Photon_Ptr->w == 0.0) {
        Photon_Ptr->dead = true;
    }
    else if (RandomNum() < CHANCE) {
        Photon_Ptr->w /= CHANCE;
    }
    else {
        Photon_Ptr->dead = true;
    }
}

/***********************************************************
 *  Fresnel reflectance.
 ****/
double RFresnel(double n1, double n2, double ca1, double* ca2_Ptr)
{
    double r;
    if (n1 == n2) {
        *ca2_Ptr = ca1; r = 0.0;
    }
    else if (ca1 > (1.0 - 1.0E-12)) {
        *ca2_Ptr = ca1;
        r = (n2 - n1) / (n2 + n1);
        r *= r;
    }
    else if (ca1 < 1.0E-6) {
        *ca2_Ptr = 0.0; r = 1.0;
    }
    else {
        double sa1 = sqrt(1 - ca1 * ca1);
        double sa2 = n1 * sa1 / n2;
        if (sa2 >= 1.0) {
            *ca2_Ptr = 0.0; r = 1.0;
        }
        else {
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

/***********************************************************
 *  Record R and T.
 ****/
void RecordR(double Refl, InputStruct* In_Ptr, PhotonStruct* Photon_Ptr, OutStruct* Out_Ptr)
{
    double x = Photon_Ptr->x, y = Photon_Ptr->y;
    short ir = (short)(sqrt(x * x + y * y) / In_Ptr->dr);
    if (ir > In_Ptr->nr - 1 || ir < 0) ir = In_Ptr->nr - 1;

    short ia = (short)(acos(-Photon_Ptr->uz) / In_Ptr->da);
    if (ia > In_Ptr->na - 1) ia = In_Ptr->na - 1;

    Out_Ptr->Rd_ra[ir][ia] += Photon_Ptr->w * (1.0 - Refl);
    Photon_Ptr->w *= Refl;
}

void RecordT(double Refl, InputStruct* In_Ptr, PhotonStruct* Photon_Ptr, OutStruct* Out_Ptr)
{
    double x = Photon_Ptr->x, y = Photon_Ptr->y;
    short ir = (short)(sqrt(x * x + y * y) / In_Ptr->dr);
    if (ir > In_Ptr->nr - 1 || ir < 0) ir = In_Ptr->nr - 1;

    short ia = (short)(acos(Photon_Ptr->uz) / In_Ptr->da);
    if (ia > In_Ptr->na - 1) ia = In_Ptr->na - 1;

    Out_Ptr->Tt_ra[ir][ia] += Photon_Ptr->w * (1.0 - Refl);
    Photon_Ptr->w *= Refl;
}

/***********************************************************
 *  Crossing decisions.
 ****/
void CrossUpOrNot(InputStruct* In_Ptr, PhotonStruct* Photon_Ptr, OutStruct* Out_Ptr)
{
    double uz = Photon_Ptr->uz;
    double uz1;
    double r = 0.0;
    short  layer = Photon_Ptr->layer;
    double ni = In_Ptr->layerspecs[layer].n;
    double nt = In_Ptr->layerspecs[layer - 1].n;

    if (-uz <= In_Ptr->layerspecs[layer].cos_crit0) r = 1.0;
    else r = RFresnel(ni, nt, -uz, &uz1);

    if (RandomNum() > r) {
        if (layer == 1) {
            Photon_Ptr->uz = -uz1;
            RecordR(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
            Photon_Ptr->dead = true;
        }
        else {
            Photon_Ptr->layer--;
            Photon_Ptr->ux *= ni / nt;
            Photon_Ptr->uy *= ni / nt;
            Photon_Ptr->uz = -uz1;
        }
    }
    else {
        Photon_Ptr->uz = -uz;
    }
}

void CrossDnOrNot(InputStruct* In_Ptr, PhotonStruct* Photon_Ptr, OutStruct* Out_Ptr)
{
    double uz = Photon_Ptr->uz;
    double uz1;
    double r = 0.0;
    short  layer = Photon_Ptr->layer;
    double ni = In_Ptr->layerspecs[layer].n;
    double nt = In_Ptr->layerspecs[layer + 1].n;

    if (uz <= In_Ptr->layerspecs[layer].cos_crit1) r = 1.0;
    else r = RFresnel(ni, nt, uz, &uz1);

    if (RandomNum() > r) {
        if (layer == In_Ptr->num_layers) {
            Photon_Ptr->uz = uz1;
            RecordT(0.0, In_Ptr, Photon_Ptr, Out_Ptr);
            Photon_Ptr->dead = true;
        }
        else {
            Photon_Ptr->layer++;
            Photon_Ptr->ux *= ni / nt;
            Photon_Ptr->uy *= ni / nt;
            Photon_Ptr->uz = uz1;
        }
    }
    else {
        Photon_Ptr->uz = -uz;
    }
}

void CrossOrNot(InputStruct* In_Ptr, PhotonStruct* Photon_Ptr, OutStruct* Out_Ptr)
{
    if (Photon_Ptr->uz < 0.0) CrossUpOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
    else                      CrossDnOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
}

/***********************************************************
 *  Move photon in glass layer.
 ****/
void HopInGlass(InputStruct* In_Ptr, PhotonStruct* Photon_Ptr, OutStruct* Out_Ptr)
{
    if (Photon_Ptr->uz == 0.0) {
        Photon_Ptr->dead = true;
    }
    else {
        StepSizeInGlass(Photon_Ptr, In_Ptr);
        Hop(Photon_Ptr);
        CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
    }
}

/***********************************************************
 *  Tissue step: hop/drop/spin with boundary handling.
 ****/
void HopDropSpinInTissue(InputStruct* In_Ptr,
    PhotonStruct* Photon_Ptr,
    OutStruct* Out_Ptr,
    int* cpt,
    double* sOut,
    short* irOut,
    short* izOut)
{
    StepSizeInTissue(Photon_Ptr, In_Ptr, Out_Ptr);

    if (HitBoundary(Photon_Ptr, In_Ptr, Out_Ptr)) {
        Hop(Photon_Ptr);
        CrossOrNot(In_Ptr, Photon_Ptr, Out_Ptr);
    }
    else {
        Hop(Photon_Ptr);
        Drop(In_Ptr, Photon_Ptr, Out_Ptr, cpt, sOut, irOut, izOut);
        Spin(In_Ptr->layerspecs[Photon_Ptr->layer].g, Photon_Ptr);
    }
}

/***********************************************************
 *  Dispatch based on layer type.
 ****/
void HopDropSpin(InputStruct* In_Ptr,
    PhotonStruct* Photon_Ptr,
    OutStruct* Out_Ptr,
    int* cpt,
    double* sOut,
    short* irOut,
    short* izOut)
{
    short layer = Photon_Ptr->layer;

    if ((In_Ptr->layerspecs[layer].mua == 0.0) && (In_Ptr->layerspecs[layer].mus == 0.0))
        HopInGlass(In_Ptr, Photon_Ptr, Out_Ptr);
    else
        HopDropSpinInTissue(In_Ptr, Photon_Ptr, Out_Ptr, cpt, sOut, irOut, izOut);

    if (Photon_Ptr->w < In_Ptr->Wth && !Photon_Ptr->dead)
        Roulette(Photon_Ptr);
}

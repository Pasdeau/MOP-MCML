/***********************************************************
 *  Numerical Recipes-like allocation helpers.
 *  Refinements: silence unused-parameter warnings.
 ****/
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

void nrerror(char error_text[])
{
    fprintf(stderr, "%s\n", error_text);
    fprintf(stderr, "...now exiting to system...\n");
    exit(1);
}

double* AllocVector(short nl, short nh)
{
    double* v;
    short i;

    v = (double*)malloc((unsigned)(nh - nl + 1) * sizeof(double));
    if (!v) nrerror("allocation failure in vector()");
    v -= nl;
    for (i = nl; i <= nh; i++) v[i] = 0.0;
    return v;
}

double** AllocMatrix(short nrl, short nrh, short ncl, short nch)
{
    short i, j;
    double** m;

    m = (double**)malloc((unsigned)(nrh - nrl + 1) * sizeof(double*));
    if (!m) nrerror("allocation failure 1 in matrix()");
    m -= nrl;

    for (i = nrl; i <= nrh; i++) {
        m[i] = (double*)malloc((unsigned)(nch - ncl + 1) * sizeof(double));
        if (!m[i]) nrerror("allocation failure 2 in matrix()");
        m[i] -= ncl;
    }

    for (i = nrl; i <= nrh; i++)
        for (j = ncl; j <= nch; j++) m[i][j] = 0.0;
    return m;
}

void FreeVector(double* v, short nl, short nh)
{
    (void)nh;
    free((char*)(v + nl));
}

void FreeMatrix(double** m, short nrl, short nrh, short ncl, short nch)
{
    (void)nch;
    short i;
    for (i = nrh; i >= nrl; i--) free((char*)(m[i] + ncl));
    free((char*)(m + nrl));
}

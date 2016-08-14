
#ifndef NANOMH_INCLUDED

#include <stdlib.h>

// data object for Gaussian proposal function
struct nMH_Gpr_ws
{
  double * sigma;
};

int nanoMH_run(size_t, int (*)(int, double *, double *, void *),
    int (*)(int, double *, double *, void *), double *, size_t,
    double *, double *, void *, void *, unsigned int);

double nMH_randn();

int nMH_Gpr(int, double *, double *, void *);

struct nMH_Gpr_ws * nMH_Gpr_make_ws(int, double *);
int nMH_Gpr_free_ws(struct nMH_Gpr_ws *);

#define NANOMH_INCLUDED
#endif


/*
  rosenbrock.c: Test nanoMH with Rosenbrock function.
  YJ Qin, AUG 16 @ Tucson
*/

#include <stdlib.h>
#include <stdio.h>
#include "nanoMH.h"

#define Sq(X) ((X) * (X))

// Rosenbrock function,
int
rosen_log_f(int N_dim, double * X, double * log_f, void * p)
{
  double x = X[0], y = X[1],
         a = ((double *)p)[0], b = ((double *)p)[1];
  double r = Sq(a - x) + b * Sq(y - Sq(x));

  * log_f = -r;

  return 0;
}

int
main(int argc, char ** argv)
{
  // Gaussian proposal function
  double sigma[2] = {0.1, 0.1};
  struct nMH_Gpr_ws * Gpr_ws_i = nMH_Gpr_make_ws(2, sigma);

  // Generate sample
  int N_pt = 30000;
  double * chain = (double *) malloc(sizeof(double) * 2 * N_pt),
         * sample = (double *) malloc(sizeof(double) * N_pt);

  double x_init[2] = {1., 1.}, rosen_par[2] = {1., 100.};
  nanoMH_run(2, & rosen_log_f, & nMH_Gpr, x_init, N_pt,
      chain, sample, (void *) rosen_par, Gpr_ws_i, 42);

  // save for plotting
  FILE * fw = fopen("chain.tmp", "wb");
  if(fw == NULL) return -1;
  fwrite(chain, sizeof(double), N_pt * 2, fw);
  fclose(fw);

  // clear and exit
  nMH_Gpr_free_ws(Gpr_ws_i);
  return 0; 
}

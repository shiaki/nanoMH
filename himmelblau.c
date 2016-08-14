
/*
  himmelblau.c: Test nanoMH with Himmelblau's function.
  YJ Qin, AUG 16 @ Tucson
*/

#include <stdlib.h>
#include <stdio.h>
#include "nanoMH.h"

#define Sq(X) ((X) * (X))

// Himmelblau's function,
int
himmel_log_f(int N_dim, double * X, double * log_f, void * p)
{
  double x = X[0], y = X[1];
  double r = Sq(Sq(x) + y - 11.) + Sq(x + Sq(y) - 7.);

  * log_f = -r;

  return 0;
}

int
main(int argc, char ** argv)
{
  // Gaussian proposal function
  double sigma[2] = {1.5, 1.5};
  struct nMH_Gpr_ws * Gpr_ws_i = nMH_Gpr_make_ws(2, sigma);

  // Generate sample
  int N_pt = 100000;
  double * chain = (double *) malloc(sizeof(double) * 2 * N_pt),
         * sample = (double *) malloc(sizeof(double) * N_pt);

  double x_init[2] = {3.5, -1.8};
  nanoMH_run(2, & himmel_log_f, & nMH_Gpr, x_init, N_pt,
      chain, sample, NULL, Gpr_ws_i, 42);

  // save for plotting
  FILE * fw = fopen("chain.tmp", "wb");
  if(fw == NULL) return -1;
  fwrite(chain, sizeof(double), N_pt * 2, fw);
  fclose(fw);
}

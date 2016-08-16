
/*
  nanoMH: simplistic implementation of the Metropolis-Hasting algorithm
  YJ Qin, AUG 2016 @ Tucson
*/

#include <stdlib.h>
#include <math.h>

#include "nanoMH.h"

#define TALLOC(type, num) ((type*) malloc(sizeof(type) * num))

int
nanoMH_run(// number of free parameters
           size_t N_dim,

           // log target function, normalization not required,
           int (* f)(int N_dim, double * x, double * log_f, void * p),

           // proposal density function, given y, generates x,
           int (* g)(int N_dim, double * x, double * y, void * p),

           // initial point,
           double * x_init,

           // size of sample,
           size_t N_pts,
           
           // chain file, size: N_size * N_dim * sizeof(double),
           double * chain,

           // density along the chain,
           double * log_f,
           
           // data object, for density function,
           void * pt_f,

           // data object, for proposal function,
           void * pt_g,
           
           // random seed,
           unsigned int seed)
{
  // status flag
  int stat;
  
  // current and proposed point, their density
  double * xc = TALLOC(double, N_dim), * xp = TALLOC(double, N_dim);
  double log_fc, log_fp;

  // N of iterations and N of sample points,
  size_t N_iter = 0, I_pt = 0;
  
  // rate of acception,
  double r_accep;

  // start: make x_init the current point, evaluate density
  for(int i = 0; i < N_dim; ++ i) xc[i] = x_init[i];
  stat = f(N_dim, xc, & log_fc, pt_f);

  // set random seed
  srand(seed);

  // start loop
  while((!stat) && (I_pt < N_pts)) 
    {
      ++ N_iter;

      // generate new point and evaluate density
      if(stat = g(N_dim, xp, xc, pt_g)) break;
      if(stat = f(N_dim, xp, & log_fp, pt_f)) break;

      // calculate acceptance ratio, and decide if accept
      if(log(1. - (double)rand() / (double)RAND_MAX) < log_fp - log_fc)
        {
          // accept this proposed point
          for(int i = 0; i < N_dim; ++ i) xc[i] = xp[i];
          for(int i = 0; i < N_dim; ++ i) chain[N_dim * I_pt + i] = xp[i];

          // write log density
          log_fc = log_fp,
          log_f[I_pt] = log_fp;

          ++ I_pt;
        }

      //printf("ITER% 8u, PTS% 8u\n\n", N_iter, I_pt);
    }

  // clean up and return stat
  free(xc), free(xp);
  return stat;
}

int
nMH_Gpr(int N_dim, double * x, double * y, void * ws)
{
  // get "kernel size" of Gaussian proposal function.
  double * sigma_i = ((struct nMH_Gpr_ws *) ws) -> sigma;

  // generate new point
  for(int i = 0; i < N_dim; ++ i) x[i] = y[i] + nMH_randn() * sigma_i[i];

  return 0;
}

struct nMH_Gpr_ws *
nMH_Gpr_make_ws(int N_dim, double * sigma)
{
  struct nMH_Gpr_ws * pt = TALLOC(struct nMH_Gpr_ws, 1);
  pt ->  sigma = TALLOC(double, N_dim);

  // copy values
  for(int i; i < N_dim; ++ i) (pt -> sigma)[i] = sigma[i];

  return pt;
}

int
nMH_Gpr_free_ws(struct nMH_Gpr_ws * pt)
{
  free(pt -> sigma), free(pt);
  return 0;
}

// Gaussian random number. Ref: Knuth, Sec. 3.4.1
double nMH_randn()
{
  static double u, v, s;
  static int p = 0;
  double x;

  if(p == 0)
    {
      do
        {
          double a = (double)rand() / RAND_MAX;
          double b = (double)rand() / RAND_MAX;

          u = 2 * a - 1;
          v = 2 * b - 1;
          s = u * u + v * v;
        }
      while(s >= 1 || s == 0);

      x = u * sqrt(-2 * log(s) / s);
    }
  else
    x = v * sqrt(-2 * log(s) / s);

  p = 1 - p;

  return x;
}

/*

  This file is a part of the wMICA R package. 
  Contact: Christoph Rau <chrau@ucla.edu>

 
  This program is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3, or (at your option)
  any later version.
 
  This program is distributed in the hope that it will be useful, but
  WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
  General Public License for more details.

*/

/* --------------------------------------------------- */

/* ---------------- vdp_mk_hp_posterior.c ------------ */

/* --------------------------------------------------- */

#include <stdio.h>
#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <R.h>
#include <Rdefines.h>
#include <time.h> 

 
/************************************************************/

/**** FUNCTIONS FOR ICMg methods ****/

void ICMgRandominit() {

  srand((unsigned)(time(0)));
  
}

/* Basic ICM with interactions only AND WEIGHTS */
void ICMgLinksIterationWeights(double *L, int *iter, int *Nlinks, int *Nnodes, int *Lindices, int *comps, int *z, double *q, int *n, double *alpha, double *beta, double *conv)
{

  /* Need to take the values from the pointer */
  int Niter = iter[0];
  int N = Nlinks[0];
  int M = Nnodes[0];
  int C = comps[0];
  double a = alpha[0];
  double b = beta[0];

  /* Auxiliary variables for computation "250"*/
  int s;
  int li;
  int l;
  int i;
  int j;
  int p;
  double w;

  /* Main iteration loop*/
  for (s=0; s<Niter; s++) {

    /* Sample new component for each link */
    for (li=0; li<N; li++) {

      l = Lindices[li]-1;
      i = L[l] -1;
      j = L[N + l] -1;
      w = L[2*N+l];


      /* Subtract the contribution of the link from the counts */
      q[i*C+z[l]-1]=q[i*C+z[l]-1]- w;
      q[j*C+z[l]-1]=q[j*C+z[l]-1]- w;
      n[z[l]-1]--;

      double uzsum = 0;
      double uz[C];

      /* Loop for computing probabilities for the components to be sampled */
      for (p=0; p<C; p++) {
	uz[p] = (n[p] + a)*(q[i*C + p]+b)*(q[j*C + p]+b)/(2*n[p]+1+M*b)/(2*n[p]+M*b);
	uzsum += uz[p];
      }

      /* Draw a new component for the links and update the counts */
      double cs = uz[0];
      int newz = -1;
      double r = rand() / ((double)(RAND_MAX) + (double)(1));
      for (p=0; p<C; p++) {

      	if (cs/uzsum >= r) {
	  newz = p;
	  conv[l] = uz[p]/uzsum;
	  break;
	}
	else
	  cs += uz[p+1];
      }

      n[newz]++;
      q[i*C+newz]=q[i*C+newz]+ w;
      q[j*C+newz]=q[j*C+newz]+ w;
      z[l] = newz+1;

    }
  }
}




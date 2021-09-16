// -*- mode: c++; c-basic-offset: 4 -*-
// (C) Nicholas Polson, James Scott, Jesse Windle, 2012-2019

// This file is part of BayesLogit.

// BayesLogit is free software: you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation, either version 3 of the License, or (at your option) any later
// version.

// BayesLogit is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

// You should have received a copy of the GNU General Public License along with
// BayesLogit.  If not, see <https://www.gnu.org/licenses/>.


#ifndef __INVERTY__
#define __INVERTY__

#include <stdio.h>

#ifdef USE_R
#include "R.h"
#endif

#include <math.h>

const double tol    = 1e-8;
const double IYPI = 3.141592653589793238462643383279502884197;

double   y_eval(double v);
void   ydy_eval(double v, double* yp, double* dyp);

double   f_eval(double v, void* params);
double  df_eval(double v, void* params);
void   fdf_eval(double v, void* params, double* fp, double* dfp);

double v_eval(double y, double tol=1e-9, int max_iter=1000);

//------------------------------------------------------------------------------

const int grid_size = 81;
const double ygrid[] = {
  0.0625,0.06698584,0.07179365,0.07694653,0.08246924,
  0.08838835,0.09473229,0.1015315,0.1088188,0.1166291,
  0.125,0.1339717,0.1435873,0.1538931,0.1649385,
  0.1767767,0.1894646,0.2030631,0.2176376,0.2332582,
  0.25,0.2679434,0.2871746,0.3077861,0.329877,
  0.3535534,0.3789291,0.4061262,0.4352753,0.4665165,
  0.5,0.5358867,0.5743492,0.6155722,0.659754,
  0.7071068,0.7578583,0.8122524,0.8705506,0.933033,
  1,1.071773,1.148698,1.231144,1.319508,
  1.414214,1.515717,1.624505,1.741101,1.866066,
  2,2.143547,2.297397,2.462289,2.639016,
  2.828427,3.031433,3.24901,3.482202,3.732132,
  4,4.287094,4.594793,4.924578,5.278032,
  5.656854,6.062866,6.498019,6.964405,7.464264,
  8,8.574188,9.189587,9.849155,10.55606,
  11.31371,12.12573,12.99604,13.92881,14.92853,
  16};

const double vgrid[] = {
  -256,-222.8609,-194.0117,-168.897,-147.0334,
  -128,-111.4305,-97.00586,-84.4485,-73.51668,
  -63.99997,-55.71516,-48.50276,-42.22387,-36.75755,
  -31.99844,-27.85472,-24.24634,-21.10349,-18.36524,
  -15.97843,-13.89663,-12.07937,-10.49137,-9.101928,
  -7.884369,-6.815582,-5.875571,-5.047078,-4.315237,
  -3.667256,-3.092143,-2.580459,-2.124095,-1.716085,
  -1.350442,-1.022007,-0.7263359,-0.4595871,-0.2184366,
  0,0.1982309,0.3784427,0.5425468,0.6922181,
  0.828928,0.953973,1.068498,1.173516,1.269928,
  1.358533,1.440046,1.515105,1.584282,1.64809,
  1.706991,1.761401,1.811697,1.858218,1.901274,
  1.941143,1.978081,2.012318,2.044068,2.073521,
  2.100856,2.126234,2.149802,2.171696,2.192042,
  2.210954,2.228537,2.244889,2.260099,2.274249,
  2.287418,2.299673,2.311082,2.321703,2.331593,
  2.340804};

//------------------------------------------------------------------------------

////////////////////////////////////////////////////////////////////////////////

//------------------------------------------------------------------------------

double y_eval(double v)
{
  double y   = 0.0;
  double r   = sqrt(fabs(v));
  if (v > tol)
    y = tan(r) / r;
  else if (v < -1*tol)
    y = tanh(r) / r;
  else
    y = 1 + (1/3) * v + (2/15) * v * v + (17/315) * v * v * v;
  return y;
}

void ydy_eval(double v, double* yp, double* dyp)
{
  // double r   = sqrt(fabs(v));

  double y = y_eval(v);
  *yp = y;

  if (fabs(v) >= tol) 
    *dyp = 0.5 * (y*y + (1-y) / v);
  else
    *dyp = 0.5 * (y*y - 1/3 - (2/15) * v);

}

double f_eval(double v, void * params)
{
  double y = *((double*) params);
  return y_eval(v) - y;
}

void fdf_eval(double v, void* params, double* fp, double* dfp)
{
  double y = *((double*)params);
  ydy_eval(v, fp, dfp);
  *fp  -= y;
}

double df_eval(double v, void * params)
{
  double f, df;
  ydy_eval(v, &f, &df);
  return df;
}

double v_eval(double y, double tol, int max_iter)
{
  double ylower = ygrid[0];
  double yupper = ygrid[grid_size-1];

  if (y < ylower) {
    return -1. / (y*y);
  } else if (y > yupper) {
    double v = atan(0.5 * y * IYPI);
    return v*v;
  }
  else if (y==1) return 0.0;
    
  double id = (log(y) / log(2.0) + 4.0) / 0.1;
  // printf("y, id, y[id], v[id]: %g, %g, %g, %g\n", y, id, ygrid[(int)id], vgrid[(int)id]);

  // C++ default is truncate decimal portion.
  int idlow  = (int)id;
  int idhigh = (int)id + 1;
  double vl  = vgrid[idlow];  // lower bound
  double vh  = vgrid[idhigh]; // upper bound

  int    iter = 0;
  double diff = tol + 1.0;
  double vnew = vl;
  double vold = vl;
  double f0, f1;

  while (diff > tol && iter < max_iter) {
    iter++;
    vold = vnew;
    fdf_eval(vold, &y, &f0, &f1);
    vnew = vold - f0 / f1;
    vnew = vnew > vh ? vh : vnew;
    vnew = vnew < vl ? vl : vnew;
    diff = fabs(vnew - vold);
    // printf("iter: %i, v: %g, diff: %g\n", iter, vnew, diff);
  }

  if (iter >= max_iter) {
    #ifndef USE_R
    fprintf(stderr, "InvertY.cpp, v_eval: reached max_iter: %i\n", iter);
    #else
    Rprintf("InvertY.cpp, v_eval: reached max_iter: %i\n", iter);
    #endif
  }

  return vnew;
}

#endif

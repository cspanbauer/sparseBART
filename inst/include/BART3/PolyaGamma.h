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

#include "simple_RNG_wrapper.h"
#include "truncated_norm.h"
#include "truncated_gamma.h"
#include "inverse_gaussian.h"

#include <cmath>
#include <vector>
#include <stdio.h>
#include <stdexcept>

using std::pow;

#ifndef __POLYAGAMMA__
#define __POLYAGAMMA__

using std::vector;

// The numerical accuracy of __PI will affect your distribution.
const double __PI = 3.141592653589793238462643383279502884197;
const double HALFPISQ = 0.5 * __PI * __PI;
const double FOURPISQ = 4 * __PI * __PI;
const double __TRUNC = 0.64;
const double __TRUNC_RECIP = 1.0 / __TRUNC;

class PolyaGamma
{

  // For sum of Gammas.
  int T;
  vector<double> bvec;

 public:

  // Constructors.
  PolyaGamma(int trunc = 200);

  // Draw.
  // double draw(double n, double z, RNG& r);
    double draw(int n, double z, rn& gen);
  double draw_sum_of_gammas(double n, double z);
    double draw_like_devroye(double z, rn& gen);

  //void draw(MF x, double a, double z, RNG& r);
  //void draw(MF x, MF a, MF z, RNG& r);

  // Utility.
  void set_trunc(int trunc);

  // Helper.
  double a(int n, double x);
  double pigauss(double x, double Z);
  double mass_texpon(double Z);
    double rtigauss(double Z, rn& gen);

  static double jj_m1(double b, double z);
  static double jj_m2(double b, double z);
  static double pg_m1(double b, double z);
  static double pg_m2(double b, double z);

};

////////////////////////////////////////////////////////////////////////////////
			       // Constructors //
////////////////////////////////////////////////////////////////////////////////

PolyaGamma::PolyaGamma(int trunc) : T(trunc), bvec(T)
{
  set_trunc(T);
} // PolyaGamma

////////////////////////////////////////////////////////////////////////////////
				 // Utility //
////////////////////////////////////////////////////////////////////////////////

void PolyaGamma::set_trunc(int trunc)
{
  
  if (trunc < 1) {
  #ifndef NTHROW
    throw std::invalid_argument("PolyaGamma(int trunc): trunc < 1.");
  #else
    #ifndef USE_R
    fprintf(stderr, "Invalid parameter: PolyaGamma(int trunc): trunc < 1.  Setting trunc=1.\n");
    #else
    Rprintf("Invalid parameter: PolyaGamma(int trunc): trunc < 1.  Setting trunc=1.\n");
    #endif
    trunc = 1;
  #endif
  }
  
  T = trunc;
  bvec.resize(T);

  for(int k=0; k < T; ++k){
    // + since we start indexing at 0.
    double d = ((double) k + 0.5);
    bvec[k] = FOURPISQ * d * d;
  }
} // set_trunc

double PolyaGamma::a(int n, double x)
{
  double K = (n + 0.5) * __PI;
  double y = 0;
  if (x > __TRUNC) {
    y = K * exp( -0.5 * K*K * x );
  }
  else if (x > 0) {
    double expnt = -1.5 * (log(0.5 * __PI)  + log(x)) + log(K) - 2.0 * (n+0.5)*(n+0.5) / x;
    y = exp(expnt);
    // y = pow(0.5 * __PI * x, -1.5) * K * exp( -2.0 * (n+0.5)*(n+0.5) / x);
    // ^- unstable for small x?
  }
  return y;
}

double PolyaGamma::pigauss(double x, double Z)
{
  double b = sqrt(1.0 / x) * (x * Z - 1);
  double a = sqrt(1.0 / x) * (x * Z + 1) * -1.0;
  double y = p_norm(b, false) + exp(2 * Z) * p_norm(a, false);
  return y;
}

double PolyaGamma::mass_texpon(double Z)
{
  double t = __TRUNC;

  double fz = 0.125 * __PI*__PI + 0.5 * Z*Z;
  double b = sqrt(1.0 / t) * (t * Z - 1);
  double a = sqrt(1.0 / t) * (t * Z + 1) * -1.0;

  double x0 = log(fz) + fz * t;
  double xb = x0 - Z + p_norm(b, true);
  double xa = x0 + Z + p_norm(a, true);

  double qdivp = 4 / __PI * ( exp(xb) + exp(xa) );

  return 1.0 / (1.0 + qdivp);
}

double PolyaGamma::rtigauss(double Z, rn& gen)
{
  Z = fabs(Z);
  double t = __TRUNC;
  double X = t + 1.0;
  if (__TRUNC_RECIP > Z) { // mu > t
    double alpha = 0.0;
    while (unif() > alpha) {
      // X = t + 1.0;
      // while (X > t)
      // 	X = 1.0 / gamma_rate(0.5, 0.5);
      // Slightly faster to use truncated normal.
      double E1 = expon_rate(1.0); double E2 = expon_rate(1.0);
      while ( E1*E1 > 2 * E2 / t) {
	E1 = expon_rate(1.0); E2 = expon_rate(1.0);
      }
      X = 1 + E1 * t;
      X = t / (X * X);
      alpha = exp(-0.5 * Z*Z * X);
    }
  }
  else {
    double mu = 1.0 / Z;
    while (X > t) {
        double Y = gen.normal(); Y *= Y;
      double half_mu = 0.5 * mu;
      double mu_Y    = mu  * Y;
      X = mu + half_mu * mu_Y - half_mu * sqrt(4 * mu_Y + mu_Y * mu_Y);
      if (unif() > mu / (mu + X))
	X = mu*mu / X;
    }
  }
  return X;
}

////////////////////////////////////////////////////////////////////////////////
				  // Sample //
////////////////////////////////////////////////////////////////////////////////

// double PolyaGamma::draw(double n, double z, RNG& r)
// {
//   return draw_sum_of_gammas(n, z, r);
// }

double PolyaGamma::draw(int n, double z, rn& gen)
{
  if (n < 1) {
  #ifndef NTHROW
    throw std::invalid_argument("PolyaGamma::draw: n < 1.");
  #else
    #ifndef USE_R
    fprintf(stderr, "PolyaGamma::draw: n < 1.  Set n = 1.\n");
    #else
    Rprintf("PolyaGamma::draw: n < 1.  Set n = 1.\n");
    #endif
    n = 1;
  #endif
  }
  double sum = 0.0;
  for (int i = 0; i < n; ++i)
      sum += draw_like_devroye(z, gen);
  return sum;
} // draw

double PolyaGamma::draw_sum_of_gammas(double n, double z)
{
  double x = 0;
  double kappa = z * z;
  for(int k=0; k < T; ++k)
    x += gamma_scale(n, 1.0) / (bvec[k] + kappa);
  return 2.0 * x;
} // draw_sum_of_gammas

double PolyaGamma::draw_like_devroye(double Z, rn& gen)
{
  // Change the parameter.
  Z = fabs(Z) * 0.5;

  // Now sample 0.25 * J^*(1, Z := Z/2).
  double fz = 0.125 * __PI*__PI + 0.5 * Z*Z;
  // ... Problems with large Z?  Try using q_over_p.
  // double p  = 0.5 * __PI * exp(-1.0 * fz * __TRUNC) / fz;
  // double q  = 2 * exp(-1.0 * Z) * pigauss(__TRUNC, Z);

  double X = 0.0;
  double S = 1.0;
  double Y = 0.0;
  // int iter = 0; If you want to keep track of iterations.

  while (true) {

    // if (unif() < p/(p+q))
    if ( unif() < mass_texpon(Z) )
      X = __TRUNC + expon_rate(1.) / fz;
    else
        X = rtigauss(Z, gen);

    S = a(0, X);
    Y = unif() * S;
    int n = 0;
    bool go = true;

    // Cap the number of iterations?
    while (go) {

      // Break infinite loop.  Put first so it always checks n==0.
      #ifdef USE_R
      if (n % 1000 == 0) R_CheckUserInterrupt();
      #endif

      ++n;
      if (n%2==1) {
	S = S - a(n, X);
	if ( Y<=S ) return 0.25 * X;
      }
      else {
	S = S + a(n, X);
	if ( Y>S ) go = false;
      }

    }
    // Need Y <= S in event that Y = S, e.g. when X = 0.

  }
} // draw_like_devroye

////////////////////////////////////////////////////////////////////////////////
			      // Static Members //
////////////////////////////////////////////////////////////////////////////////

double PolyaGamma::jj_m1(double b, double z) 
{
    z = fabs(z);
    double m1 = 0.0;
    if (z > 1e-12)
	m1 = b * tanh(z) / z;
    else
	m1 = b * (1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6));
    return m1;
}

double PolyaGamma::jj_m2(double b, double z)
{
    z = fabs(z);
    double m2 = 0.0;
    if (z > 1e-12)
	m2 = (b+1) * b * pow(tanh(z)/z,2) + b * ((tanh(z)-z)/pow(z,3));
    else
	m2 = (b+1) * b * pow(1 - (1.0/3) * pow(z,2) + (2.0/15) * pow(z,4) - (17.0/315) * pow(z,6), 2) +
	    b * ((-1.0/3) + (2.0/15) * pow(z,2) - (17.0/315) * pow(z,4));
    return m2;
}

double PolyaGamma::pg_m1(double b, double z)
{
    return jj_m1(b, 0.5 * z) * 0.25;
}
 
double PolyaGamma::pg_m2(double b, double z)
{
    return jj_m2(b, 0.5 * z) * 0.0625;
}

#endif
////////////////////////////////////////////////////////////////////////////////
				 // APPENDIX //
////////////////////////////////////////////////////////////////////////////////

// It was not faster to use "vectorized" versions of r.gamma or r.igamma.

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



// MAKE SURE YOU CALL GetRNGSeed() and PutRNGSeed() WHEN USING THESE FUNCTIONS!!!

#include "R.h"
#include "Rmath.h"
#include <cmath>

// Define MACROS

#ifndef SIMPLE_RNG_WRAPPER
#define SIMPLE_RNG_WRAPPER

#define RCHECK 1000

const double SQRT2PI = 2.50662827;

inline void check_R_interupt(int count)
{
    #ifdef USE_R
    if (count % RCHECK == 0) R_CheckUserInterrupt();
    #endif
}

// CDF
#define p_norm(X, USE_LOG) R::pnorm((X), 0.0, 1.0, true, (USE_LOG))
#define p_gamma_scale(X, SHAPE, SCALE, USE_LOG) R::pgamma((X), (SHAPE), (SCALE), 1, (USE_LOG))
#define p_gamma_rate(X, SHAPE, RATE, USE_LOG) R::pgamma((X), (SHAPE), (1./(RATE)), 1, (USE_LOG))

// Random variates
#define expon_mean(MEAN) R::rexp((MEAN))
#define expon_rate(RATE) R::rexp((1./(RATE)))
#define unif() R::runif(0.0, 1.0)
#define norml(MEAN, SD) R::rnorm((MEAN), (SD))
#define gamma_scale(SHAPE, SCALE) R::rgamma((SHAPE), (SCALE))
#define gamma_rate(SHAPE, RATE) R::rgamma((SHAPE), (1./(RATE)))

#endif
// Output



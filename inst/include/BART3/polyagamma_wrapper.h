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



// MAKE SURE YOU USE GetRNGState() and PutRNGState() !!!

#include "PolyaGamma.h"
#include "PolyaGammaApproxSP.h"
#include "simple_RNG_wrapper.h"

#ifndef POLYAGAMMA_HYBRID
#define POLYAGAMMA_HYBRID

extern void rpg_hybrid(double& x, double h, double z, rn& gen)
{
    PolyaGamma dv(1000);
    // PolyaGammaApproxAlt alt;
    PolyaGammaApproxSP sp;

#ifdef USE_R
    GetRNGstate();
#endif
    if (h > 170) {
        double m = dv.pg_m1(h,z);
        double v = dv.pg_m2(h,z) - m*m;
        x = m + sqrt(v)*gen.normal();
    }
    else if (h > 13) {
        sp.draw(x, h, z, gen);
    }
    else if (h==1 || h==2) {
        x = dv.draw((int)h, z, gen);
    }
	// Need to review "alt" sampler.
    // else if (b > 1) {
    //     x[i] = alt.draw(b, z[i]);
    // }
    else if (h > 0) {
        x = dv.draw_sum_of_gammas(h, z);
    }
    else {
        x = 0.0;
    }

#ifdef USE_R
    PutRNGstate();
#endif
}

#endif

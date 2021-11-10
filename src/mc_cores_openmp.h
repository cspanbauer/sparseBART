/*
 *  sparseBART: sparse Bayesian Additive Regression Trees
 *  Copyright (C) 2021 Charles Spanbauer
 *
 *  This file is part of sparseBART.
 *
 *  sparseBART is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; version 3 of the License, or
 *  (at your option) any later version.
 *
 *  sparseBART is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with sparseBART; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-3
 */

#ifndef NoRcpp

RcppExport SEXP mc_cores_openmp() {

#else

int mc_cores_openmp() {

#endif

#ifdef _OPENMP

int mc_cores_openmp=omp_get_num_threads();

#else

int mc_cores_openmp=0;

#endif

#ifndef NoRcpp

return Rcpp::wrap(mc_cores_openmp);

#else

return mc_cores_openmp;

#endif

}


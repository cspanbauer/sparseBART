/*
 *  BART: Bayesian Additive Regression Trees
 *  Copyright (C) 2017-2018 Robert McCulloch, Rodney Sparapani
 *                          and Charles Spanbauer
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, a copy is available at
 *  https://www.R-project.org/Licenses/GPL-2
 */


#ifndef GUARD_vs_h
#define GUARD_vs_h

#include "RcppEigen.h"

typedef Eigen::Map<Eigen::MatrixXd> mm;
typedef Eigen::VectorXd vm;


void unscented(vm& _xi_prior_mean, mm& _xi_prior_var, vm& _psi_prior_mean, vm& _xi_prior_var){
  size_t pm1=_psi_prior_mean.size()-1;
  Eigen::EigenSolver<MatrixXd> es;
  es.compute(_psi_prior_mean, true);
  double alpha=0.001; double lambda=pm1*(alpha*alpha-1);
  vm wts (2*pm1+1,0.);
  mm pts;
  wts[0] = lambda/(pm1*alpha*alpha);
  pts[0] = _psi_prior_mean;
  _xi_prior_mean.setZero();
  _xi_prior_mean = wts[0]*pts[0];
  _xi_prior_var.setZero();
  for(size_t j=1;j<wts.size();j++){
    wts[j] = 1/(2*alpha*alpha);
  }
  for(size_t j=1;j<(2*pm1+1);j++){
    pts[j] = _psi_prior_mean+pow(pm1+lambda,0.5)*es.eigenvalues()[j-1]*es.eigenvectors().col(j-1);
  }
  for(size_t j=(2*pm+1);j<wts.size();j++){
    pts[j] = _psi_prior_mean-pow(pm1+lambda,0.5)*es.eigenvalues()[j-pm1-1]*es.eigenvectors().col(i-pm1-1);
  }
  vm z1 (2*pm1+1,0.);
  vm z2 (
  
  _xi_prior_var = (1-alpha*alpha+2)*(
}

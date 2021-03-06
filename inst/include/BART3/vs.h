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


#ifndef GUARD_vs_h
#define GUARD_vs_h

#include "polyagamma_wrapper.h"
#include "inverse_gaussian.h"
#include "RcppEigen.h"

typedef Eigen::MatrixXd mm;
typedef Eigen::VectorXd vm;

class vs {
public:
  vs(size_t _sparse, size_t _p);
  ~vs();
  vs& operator=(const vs&);
  void turn_on();
  bool is_on();
  size_t get_sparse();
  void setvs_dart(double _a, double _b, int *_grp, bool _aug,
                  double _rho=0., double _theta=0., double _omega=1.);
  void setvs_dart(double _a, double _b, bool _aug, double _rho, double _theta, double _omega);
  void setvs_ln(double _bprvar, double _tprprm, mm _an, size_t _noan, bool _lasso, rn& gen);
  void draw_s_dart(std::vector<double>& lpv, std::vector<size_t>& nv, rn& gen);
  void draw_s_grp_dart(std::vector<size_t>& nv, std::vector<double>& lpv, rn& gen, double * grp);
  void draw_theta0_dart(std::vector<double>& lpv,rn& gen);
  void draw_s_ln(std::vector<size_t>& nv, std::vector<double>& lpv, rn &gen);
  void draw_tau2_ln(rn &gen);
  void draw_lambda_lnLASSO(rn& gen);
  void set_R(size_t _R);
protected:
  size_t sparse;
  size_t p;
  bool isOn, const_theta;
  // DART parameters
  double a,b,rho;
  bool aug;
  double theta, omega;
  // Logit-Normal parameters
  mm an;
  mm invD;
  double bprvar,tprprm;
  size_t noan;
  vm beta;
  vm psi;
  vm sig2Inv;
  vm resid;
  double lse;
  double tau2;
  double *grp;
  size_t R;
  double max_psi;
  bool lasso;
  double lambda;
};

vs::vs(size_t _sparse, size_t _p):grp(0)
{
  this->isOn=false;
  this->sparse=_sparse;
  this->p=_p;
  cout << "P: " << p << '\n';
}

vs::~vs()
{
  if(grp) delete[] grp;
}
  
vs& vs::operator=(const vs&)
  {
    return *this;
  }
void vs::turn_on() {
  this->isOn=true;
}
bool vs::is_on()
{
  return(isOn);
}

size_t vs::get_sparse()
{
  return sparse;
}

void vs::setvs_dart(double _a, double _b, int *_grp, bool _aug,
                    double _rho, double _theta, double _omega)
{
 grp = new double[p];
 for(size_t i=0;i<p;++i) grp[i]=_grp[i];
 if(_rho==0.) for(size_t i=0;i<p;++i) _rho += 1./grp[i];
 this->setvs_dart(_a, _b, _aug, _rho, _theta, _omega);
}
void vs::setvs_dart(double _a, double _b, bool _aug, double _rho, double _theta, double _omega)
{
 this->a=_a; this->b=_b; this->rho=_rho; this->aug=_aug; this->omega=_omega;
 if(_theta==0.){
   this->const_theta=false;
   this->theta=1.;
 }
 else {
   this->const_theta=true;
   this->theta=_theta;
 }
}
void vs::setvs_ln(double _bprvar, double _tprprm, mm _an, size_t _noan, bool _lasso, rn& gen)
{
 this->bprvar=_bprvar; this->tprprm=_tprprm;
 this->an=_an; this->noan=_noan; this->lasso=_lasso;
 beta.resize(noan);
 sig2Inv.resize(noan);
 resid.resize(p);
 double sum_absBeta=0.;
 for(size_t t=0;t<noan;t++) {
   beta[t]=0.01*gen.normal();
   sig2Inv[t]=1/(beta[t]*beta[t]);
   sum_absBeta+=abs(beta[t]);
 }
 psi.resize(p);
 double sum_resid2=0.;
 for(size_t j=0;j<p;j++) {
   psi[j]=0.;
   resid[j]=psi[j]-an.row(j)*beta;
   sum_resid2+=pow(resid[j],2.);
 }
 this->tau2=1.;
 this->lse=::log(p);
 this->max_psi=psi[0];
 this->lambda=(double)(p-1)*pow(sum_resid2/(double)(p-1),0.5)/sum_absBeta;
}

//--------------------------------------------------
//draw variable splitting probabilities from Dirichlet (Linero, 2018)
void vs::draw_s_dart(std::vector<double>& lpv, std::vector<size_t>& nv, rn& gen){
// Now draw s, the vector of splitting probabilities
  std::vector<double> _theta(p);
  for(size_t j=0;j<p;j++) _theta[j]=theta/(double)p+(double)nv[j];
  //gen.set_alpha(_theta);
  lpv=gen.log_dirichlet(_theta);
}

void vs::draw_s_grp_dart(std::vector<size_t>& nv, std::vector<double>& lpv, rn& gen, double *grp){
  size_t p=nv.size();
// Now draw s, the vector of splitting probabilities
  std::vector<double> _theta(p);
  for(size_t j=0;j<p;j++) {
    if(grp) _theta[j]=theta/(rho*grp[j])+(double)nv[j];
    else _theta[j]=theta/rho+(double)nv[j];
  }
  //gen.set_alpha(_theta);
  lpv=gen.log_dirichlet(_theta);
}

void vs::draw_theta0_dart(std::vector<double>& lpv,rn& gen){
  // Draw sparsity parameter theta_0 (Linero calls it alpha); see Linero, 2018
  // theta / (theta + rho ) ~ Beta(a,b)
  // Set (a=0.5, b=1) for sparsity
  // Set (a=1, b=1) for non-sparsity
  // rho = p usually, but making rho < p increases sparsity
  if(!const_theta){
    double sumlpv=0.,lse;
    
    std::vector<double> lambda_g (1000,0.);
    std::vector<double> theta_g (1000,0.);
    std::vector<double> lwt_g (1000,0.);
    for(size_t j=0;j<p;j++) sumlpv+=lpv[j];
    for(size_t k=0;k<1000;k++){
      lambda_g[k]=(double)(k+1)/1001.;
      theta_g[k]=(lambda_g[k]*rho)/(1.-lambda_g[k]);
      double theta_log_lik=lgamma(theta_g[k])-(double)p*lgamma(theta_g[k]/(double)p)+(theta_g[k]/(double)p)*sumlpv;
      double beta_log_prior=(a-1.)*log(lambda_g[k])+(b-1.)*log(1.-lambda_g[k]);
//      cout << "SLP: " << sumlogpv << "\nTLL: " << theta_log_lik << "\nBLP: " << beta_log_prior << '\n';
      lwt_g[k]=theta_log_lik+beta_log_prior;      
    }
    lse=log_sum_exp(lwt_g);
    for(size_t k=0;k<1000;k++) {
      lwt_g[k]=exp(lwt_g[k]-lse);
//      cout << "LWT: " << lwt_g[k] << '\n';
    }
    gen.set_wts(lwt_g);    
    theta=theta_g[gen.discrete()];
  } 
}

void vs::draw_s_ln(std::vector<size_t>& nv, std::vector<double>& lpv, rn &gen)
{
  // Draw psi
  double psi_prior_mean=0.;
  for(size_t j=0;j<(p-1);j++){
    //    this->lse=max_psi+::log(::exp(lse)-::exp(psi[j]-max_psi));
    this->lse=::log(::exp(lse)-::exp(psi[j]));
    double xi=psi[j]-lse;
    double kappa=(double)nv[j]-0.5*(double)R;
    double omega;
    // Draw Polya-gamma
    rpg_hybrid(omega,(double)R,xi,gen);
    // Draw psi[j]
    if(noan!=0) psi_prior_mean = an.row(j)*beta;
    double psiVar = tau2/(1.+tau2*omega);
    double psiMean = (psi_prior_mean+tau2*(kappa+omega*lse))/(1+tau2*omega);
    psi[j]=pow(psiVar,0.5)*gen.normal()+psiMean;
    if(psi[j]>max_psi) max_psi=psi[j];
    //    this->lse=max_psi+::log(::exp(this->lse)+::exp(psi[j]-max_psi));
    this->lse=::log(::exp(this->lse)+::exp(psi[j]));
  }
  // Draw annotation beta
  if(noan!=0&!lasso){
    mm betaVar = ((bprvar*an.transpose()*an+tau2*Eigen::MatrixXd::Identity(noan,noan))/(tau2*bprvar)).inverse();
    vm betaMean = betaVar*((1.0/tau2)*an.transpose()*psi);
    vm Zscores;
    Zscores.resize(noan);
    for(size_t t=0;t<noan;t++){
      Zscores[t]=gen.normal();
    }
    Eigen::LLT<Eigen::MatrixXd> LLTbetaVar(betaVar);
    this->beta = LLTbetaVar.matrixL()*Zscores+betaMean;
  }
  else if(noan!=0&lasso){
    this->invD = sig2Inv*Eigen::MatrixXd::Identity(noan,noan);
    mm invA = (an.transpose()*an+invD).inverse();
    vm betaMean = invA*an.transpose()*psi;
    mm betaVar = (resid.transpose()*resid)/(double)(p-1);
    vm Zscores;
    Zscores.resize(noan);
    for(size_t t=0;t<noan;t++){
      Zscores[t]=gen.normal();
    } 
    Eigen::LLT<Eigen::MatrixXd> LLTbetaVar(betaVar);
    this->beta = LLTbetaVar.matrixL()*Zscores+betaMean;
    for(size_t j=0;j<p-1;j++){
      resid[j]=psi[j]-an.row(j)*beta;
    }
  }
  // Compute log(s[j]) from psi
  for(size_t j=0;j<p;j++){
    lpv[j]=psi[j]-lse;
  }
}

void vs::draw_tau2_ln(rn& gen)
{
  double sum_eta2=0.;
  if(noan==0){
    for(size_t j=0;j<p;j++){
      sum_eta2+=pow(psi[j],2.);
    }
    double phi = 1.0/gen.gamma(2.,3./tau2+1./(tprprm*tprprm));
    this->tau2=1.0/gen.gamma(0.5*double(p+3),0.5*sum_eta2+3./phi);
  }
  if(noan!=0&!lasso) {
    for(size_t j=0;j<p;j++){
      sum_eta2+=pow(psi[j]-an.row(j)*beta,2.);
    }
    double phi = 1.0/gen.gamma(2.,3./tau2+1./(tprprm*tprprm));
    this->tau2=1.0/gen.gamma(0.5*double(p+3),0.5*sum_eta2+3./phi);
  }
  else if(noan!=0&lasso) {
    for(size_t j=0;j<p;j++){
      sum_eta2+=pow(resid[j],2.);
    }
    this->tau2=1.0/gen.gamma(0.5*double(p-2+noan),2.0/(sum_eta2+beta.transpose()*invD*beta));
  }
}

void vs::draw_lambda_lnLASSO(rn& gen) {
  double sum_oneOversig2Inv=0.;
  for(size_t t=0;t<noan;t++){
    sig2Inv[t]=igauss(pow(lambda*lambda*tau2/(beta[t]*beta[t]),0.5),lambda*lambda,gen);
    sum_oneOversig2Inv+=1.0/sig2Inv[t];
  }
  double shp=1.+0.5*(double)(noan);
  double scl=0.1+0.5*sum_oneOversig2Inv;
  this->lambda=gen.gamma(shp,1.0/scl);
}

void vs::set_R(size_t _R)
{
  this->R=_R;
}

#endif

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

#ifndef GUARD_heterbart_h
#define GUARD_heterbart_h

#include "bart.h"
#include "heterbartfuns.h"
#include "heterbd.h"

class heterbart : public bart
{
  public:
   heterbart():bart() { }
   heterbart(size_t m):bart(m) { }
   void pr();
   void draw(double *sigma, rn& gen, int shards=1);
};

//--------------------------------------------------
void heterbart::pr()
{
   cout << "+++++heterbart object:\n";
   bart::pr();
}
//--------------------------------------------------
void heterbart::draw(double *sigma, rn& gen, int shards)
{
   size_t i=0;
   for(size_t j=0;j<m;j++) {
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) {
         allfit[k] = allfit[k]-ftemp[k];
         r[k] = y[k]-allfit[k];
      }
      if(heterbd(t[j],xi,di,pi,sigma,nv,pv,false,gen,shards)) i++;
      heterdrmu(t[j],xi,di,pi,sigma,gen);
      fit(t[j],xi,p,n,x,ftemp);
      for(size_t k=0;k<n;k++) allfit[k] += ftemp[k];
   }
   //   accept=i/(double)m;
}

#endif

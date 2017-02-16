/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:17:49
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 23:46:23
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./projector.h"

namespace var {

void WavefunProjector::init(const input::Parameters& parms) 
{
  int info;
  gutzwiller_proj_ = parms.set_value("gutzwiller_proj", false, info);
  double g = 1.0;
  if (gutzwiller_proj_) g = parms.set_value("g", 1.0);
  pfactors_.insert({"g", g});
  // gw ratio
  gw_ratio_.resize(3);
  gw_ratio_[0] = 1.0/g;  // nd_increament = -1
  gw_ratio_[1] = 1.0;    // nd_increament = 0
  gw_ratio_[2] = g;      // nd_increament = 1
}

} // end namespace var

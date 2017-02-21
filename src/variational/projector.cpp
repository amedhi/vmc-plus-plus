/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:17:49
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-21 10:08:14
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./projector.h"

namespace var {

void WavefunProjector::init(const input::Parameters& parms) 
{
  int info;
  gutzwiller_proj_ = parms.set_value("gutzwiller_proj", false, info);
  double g = 1.0;
  if (gutzwiller_proj_) {
    g = parms.set_value("gfactor", 1.0);
    varparms_.add("g", g, 1.0E-6, 1.0);
  }
  //pfactors_.insert({"g", g});
  // gw ratio
  gw_ratio_ = {1.0,1.0,1.0};
  //gw_ratio_[0] = 1.0/g;  // nd_increament = -1
  //gw_ratio_[1] = 1.0;    // nd_increament = 0
  //gw_ratio_[2] = g;      // nd_increament = 1
}

void WavefunProjector::update(const input::Parameters& inputs) 
{ 
  varparms_.update(inputs); 
  if (gutzwiller_proj_) set_gw_ratio();
} 

void WavefunProjector::update(const std::vector<double>& vparms, const unsigned& begin,
    const unsigned& end) 
{ 
  varparms_.update(vparms,begin,end); 
  if (gutzwiller_proj_) set_gw_ratio();
}

void WavefunProjector::set_gw_ratio(void) 
{ 
  int pos = varparms_.at("g");
  double g = varparms_.values()[pos];
  //std::cout << "### g = " << g << "\n";
  gw_ratio_[0] = 1.0/g;  // nd_increament = -1
  gw_ratio_[1] = 1.0;    // nd_increament = 0
  gw_ratio_[2] = g;      // nd_increament = 1
} 




} // end namespace var

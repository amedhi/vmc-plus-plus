/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:03:44
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 23:46:15
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef PROJECTOR_H
#define PROJECTOR_H

#include <iostream>
#include <vector>
#include <map>
#include "../scheduler/worker.h"
//#include "../lattice/lattice.h"
//#include "../lattice/graph.h"
//#include "../variational/matrix.h"

namespace var {

class WavefunProjector 
{
public:
  WavefunProjector() {}
  WavefunProjector(const input::Parameters& parms) { init(parms); }
  ~WavefunProjector() {}
  void init(const input::Parameters& parms); 
  const double& gw_ratio(const int& nd_incre) const { return gw_ratio_[nd_incre+1]; } 
private:
  using name_value_pair = std::pair<std::string,double>;
  bool gutzwiller_proj_{false};
  //int num_gw_factors_{0};
  std::vector<double> gw_ratio_;
  std::map<std::string,double> pfactors_;
};


} // end namespace var

#endif
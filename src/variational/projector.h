/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:03:44
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-25 11:45:47
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef PROJECTOR_H
#define PROJECTOR_H

#include <iostream>
#include <vector>
#include <map>
#include "../scheduler/worker.h"
#include "./varparm.h"
//#include "../lattice/lattice.h"
//#include "../lattice/graph.h"
//#include "../variational/matrix.h"

namespace var {

enum class pp {gutzwiller, end};

class WavefunProjector 
{
public:
  WavefunProjector() {}
  WavefunProjector(const input::Parameters& parms) { init(parms); }
  ~WavefunProjector() {}
  void init(const input::Parameters& inputs); 
  void update(const input::Parameters& inputs); 
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0);
  const bool& have_gutzwiller(void) const { return gutzwiller_proj_; }
  double gw_factor(void) const; 
  const double& gw_ratio(const int& nd_incre) const { return gw_ratio_[nd_incre+1]; } 
  const VariationalParms& varparms(void) const { return varparms_; }
private:
  using vparm_t = std::pair<std::string,double>;
  bool gutzwiller_proj_{false};
  //int num_gw_factors_{0};
  std::vector<double> gw_ratio_;
  //std::map<std::string,double> pfactors_;
  VariationalParms varparms_;

  void set_gw_ratio(void);
};


} // end namespace var

#endif
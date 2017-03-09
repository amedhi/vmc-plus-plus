/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:07:37
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-09 17:01:24
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef STOCHASTIC_RECONF_H
#define STOCHASTIC_RECONF_H

#include <iostream>
#include "../utils/utils.h"
#include "./simulator.h"

namespace vmc {

class StochasticReconf 
{
public:
  StochasticReconf() {} 
  StochasticReconf(const input::Parameters& parms); 
  ~StochasticReconf() {}
  int init(const input::Parameters& parms, const Simulator& simulator);
  int optimize(Simulator& simulator);
  //const var::parm_vector& vp(void) { return varparms; }
private:
  Observable optimal_parms_;
  unsigned num_parms_;
  var::parm_vector vparms_;
  var::parm_vector lbound_;
  var::parm_vector ubound_;
  Eigen::MatrixXd sr_matrix_;
  Eigen::VectorXd grad_;
  // Mann-Kendall trend test for converegence
  util::MK_Statistic mk_statistic_;

  // optimization parameters
  unsigned num_sim_samples_{1000};
  unsigned num_opt_samples_{30};
  unsigned max_iter_{500};
  unsigned refinement_cycle_{100};
  unsigned mk_series_len_{40};
  double start_tstep_{0.05};
  double mk_thresold_{0.30};
};


} // end namespace vmc

#endif
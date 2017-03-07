/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef VMC_H
#define VMC_H

#include <iostream>
#include "./simulator.h"

namespace vmc {

class VMC : public scheduler::Worker
{
public:
  VMC(const input::Parameters& parms); 
  ~VMC() {}
  int start(input::Parameters& parms) override { return 0; }
  int run(input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
  const Simulator& sim(void) { return simulator; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  Simulator simulator;
  //var::parm_vector varparms;
  var::parm_vector varp_lb_;
  var::parm_vector varp_ub_;

  // optimization parameters
  bool optimization_mode_{false};
  unsigned num_opt_samples_{1};
  Observable opt_varp_;

  // stochastic reconfiguration parameters
  unsigned sr_max_iter_{500};
  double sr_tstep_{0.05};
  unsigned sr_max_mklen_{30};
  double sr_mktrend_tol_{0.30};

  int run_optimization(input::Parameters& inputs);
  static double enfunc(const std::vector<double>& x, std::vector<double>& grad, 
    void *my_func_data);
};

double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data);

} // end namespace monte carlo

#endif

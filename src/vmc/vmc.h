/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef VMC_H
#define VMC_H

#include <iostream>
#include "./simulator.h"
#include "./stochastic_reconf.h"

namespace vmc {

class VMC : public scheduler::Worker
{
public:
  VMC(const input::Parameters& parms); 
  ~VMC() {}
  int start(const input::Parameters& parms) override { return 0; }
  int run(const input::Parameters& parms) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
  const Simulator& sim(void) { return simulator; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  Simulator simulator;
  StochasticReconf sreconf;
  bool optimization_mode_{false};

  //static double enfunc(const std::vector<double>& x, std::vector<double>& grad, 
  //  void *my_func_data);
};

//double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data);

} // end namespace vmc

#endif

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
  int start(input::Parameters& parms) override;
  void run(void) override {} 
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
private:
  Simulator simulator;
  var::parm_vector varparms;
  var::parm_vector vp_lbound;
  var::parm_vector vp_ubound;
};


} // end namespace monte carlo

#endif

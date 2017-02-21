/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-21 09:49:53
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "vmc.h"

namespace vmc {

VMC::VMC(const input::Parameters& parms) : simulator(parms)
{
}

int VMC::start(input::Parameters& parms) 
{
  std::cout << " Simulator::run\n";
  //simulator.run(parms);
  simulator.get_variational_parms(varparms_);
  std::cout << "var parms = " << varparms_.size();
  int j = 0;
  for (int i=0; i<10; ++i) {
    simulator.run(varparms_);
    varparms_[j++] += 0.2;
    if (j == varparms_.size()) j = 0;
  }

  return 0;
}

void VMC::print_copyright(std::ostream& os)
{
  Simulator::print_copyright(os);
}


} // end namespace mc

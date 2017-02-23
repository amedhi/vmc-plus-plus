/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-24 00:06:15
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "vmc.h"

namespace vmc {

VMC::VMC(const input::Parameters& inputs) : simulator(inputs)
{
}

int VMC::start(input::Parameters& inputs) 
{
  // starting
  simulator.init(inputs);

  // normal run
  if (!simulator.optimizing_mode()) {
    std::cout << " starting vmc run\n";
    simulator.run();
    simulator.print_results();
    return 0;
  }

  // optimizing run
  std::cout << " starting vmc optimization\n";
  simulator.get_variational_parms(varparms);
  //simulator.run(varparms, true);
  return 0;
  /*
  std::cout << "var parms = " << varparms_.size();
  simulator.optimizing(varparms_,energy,energy_grad);
  int j = 0;
  for (int i=0; i<10; ++i) {
    simulator.run(varparms_,energy,energy_grad);
    varparms_[j++] += 0.2;
    if (j == varparms_.size()) j = 0;
  }
  return 0;
  */
}

void VMC::print_copyright(std::ostream& os)
{
  Simulator::copyright_msg(os);
}


} // end namespace mc

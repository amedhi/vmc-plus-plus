/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-18 23:14:05
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
  simulator.run(parms);

  return 0;
}

void VMC::print_copyright(std::ostream& os)
{
  Simulator::print_copyright(os);
  //os << "#" << std::string(72,'-') << "\n";
  //os << "#" << " Program: VMC Simulation\n";
  //os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  //os << "#" << std::string(72,'-') << "\n";
}


} // end namespace mc

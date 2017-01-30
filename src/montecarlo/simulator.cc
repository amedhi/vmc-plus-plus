/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   amedhi
* Last Modified time: 2017-01-30 20:45:54
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "simulator.h"

namespace mc {

VMC::VMC(input::Parameters& parms) 
  : LatticeGraph(parms) 
  , Model(parms, LatticeGraph::lattice()) 
  , psi_var_(parms, LatticeGraph::lattice())
{
  //Model::construct(LatticeGraph::lattice(), parms);
}

int VMC::start(input::Parameters& parms) 
{
  // this function must be override
  std::cout << "Simulator::start\n";
  //lattice::Lattice lattice(parms);
  std::cout << "lattice: name = " << lattice().name() << "\n";
  std::cout << "lattice: num sites = " << lattice().num_sites() << "\n";

  return 0;
}

void VMC::print_copyright(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}


} // end namespace mc

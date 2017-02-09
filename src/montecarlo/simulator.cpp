/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-01 22:20:53
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "simulator.h"

namespace mc {

VMC::VMC(input::Parameters& parms) 
  : graph_(parms) 
  , model_(parms, graph_.lattice()) 
  , wavefunc_(parms, graph_)
{
  //std::cout << "-----------hi----------\n";
  //Model::construct(LatticeGraph::lattice(), parms);
}

int VMC::start(input::Parameters& parms) 
{
  // this function must be override
  std::cout << "Simulator::start\n";
  //lattice::Lattice lattice(parms);
  std::cout << "lattice: name = " << graph_.lattice().name() << "\n";
  std::cout << "lattice: num sites = " << graph_.lattice().num_sites() << "\n";

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

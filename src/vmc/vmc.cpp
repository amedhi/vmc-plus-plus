/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-02 23:57:29
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <nlopt.hpp>
#include "vmc.h"
#include "../optimizer/LBFGS.h"

namespace vmc {

VMC::VMC(const input::Parameters& inputs) : simulator(inputs)
{
}

int VMC::start(input::Parameters& inputs) 
{
  // starting
  simulator.start(inputs);

  // normal run
  if (!simulator.optimizing_mode()) {
    std::cout << " starting vmc run\n";
    simulator.run();
    simulator.print_results();
    return 0;
  }

  // optimizing run
  std::cout << " starting vmc optimization\n";
  simulator.get_vparm_values(varparms);

  nlopt::opt opt(nlopt::LD_MMA, varparms.size());


  /*using namespace LBFGSpp;
  LBFGSParam<double> lbfgs_param;
  LBFGSSolver<double> solver(lbfgs_param);
  //simulator.run(varparms, true);
  double emin;
  int niter = solver.minimize(simulator, varparms, emin);
  std::cout << niter << " iterations" << std::endl;
  std::cout << "x = \n" << varparms.transpose() << std::endl;
  std::cout << "f(x) = " << emin << std::endl;
  */
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
  return 0;
}

void VMC::print_copyright(std::ostream& os)
{
  Simulator::copyright_msg(os);
}


} // end namespace mc

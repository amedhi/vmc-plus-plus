/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-10 15:42:21
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "vmc.h"
//#include <nlopt.hpp>
//#include "../optimizer/LBFGS.h"

namespace vmc {

VMC::VMC(const input::Parameters& inputs) : simulator(inputs)
{
  optimization_mode_ = inputs.set_value("optimizing_run",false);
  if (optimization_mode_) sreconf.init(inputs, simulator);
}

int VMC::run(const input::Parameters& inputs) 
{
  // optimization run
  if (optimization_mode_) {
    simulator.start(inputs, true, true);
    if (sreconf.optimize(simulator)) {
      simulator.run_opt_simulation(sreconf.optimal_parms());
      simulator.print_results();
    }
    return 0;
  }

  // normal run
  std::cout << " starting vmc run\n";
  simulator.start(inputs);
  simulator.run_simulation();
  simulator.print_results();
  return 0;
}

  /*
  //nlopt::opt opt(nlopt::LD_MMA, varparms.size());
  nlopt::opt opt(nlopt::LD_LBFGS, varparms.size());
  opt.set_lower_bounds(simulator.vparm_lbound());
  opt.set_upper_bounds(simulator.vparm_ubound());
  opt.set_xtol_rel(1e-3);
  opt.set_min_objective(wrapper, &simulator);
  double mine;
  nlopt::result result = opt.optimize(varparms, mine);
  //std::cout << nlopt::result << "\n";
  */

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

/*
double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data)
{
  Eigen::VectorXd en_grad(x.size());
  double en = reinterpret_cast<Simulator*>(my_data)->energy_function(x,en_grad);
  for (int i=0; i<grad.size(); ++i) grad[i] = en_grad[i];
  return en;
}
*/

void VMC::print_copyright(std::ostream& os)
{
  Simulator::copyright_msg(os);
}


} // end namespace mc

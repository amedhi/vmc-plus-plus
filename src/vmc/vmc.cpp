/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-08 00:05:08
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <nlopt.hpp>
#include "vmc.h"
#include "../utils/utils.h"
//#include "../optimizer/LBFGS.h"

namespace vmc {

VMC::VMC(const input::Parameters& inputs) : simulator(inputs)
{
  optimization_mode_ = inputs.set_value("optimizing_run",false);
}

int VMC::run(input::Parameters& inputs) 
{
  // optimization run
  if (optimization_mode_) return run_optimization(inputs);
  // normal run
  std::cout << " starting vmc run\n";
  simulator.start(inputs);
  simulator.run();
  simulator.print_results();
  return 0;
}

int VMC::run_optimization(input::Parameters& inputs)
{
  // varp observable 
  //opt_varp_.init(inputs, "OptParams", simulator.varp_names());
  //opt_varp_.set_elements(simulator.varp_names());
  //opt_varp_.switch_on();
  // optimizing run
  std::cout << " starting vmc optimization\n";
  simulator.start(inputs, true);
  // variational parameters bound
  //varparms = simulator.varp_values();
  varp_lb_ = simulator.varp_lbound();
  varp_ub_ = simulator.varp_ubound();
  var::parm_vector varp(varp_lb_.size());
  // Stochastic Reconfiguration matrix
  Eigen::MatrixXd sr_matrix;
  Eigen::VectorXd grad;
  // Mann-Kendall statistic
  util::MK_Statistic mk_statistic(varp.size(), sr_max_mklen_);
  // start optimization
  for (unsigned n=0; n<num_opt_samples_; ++n) {
    std::cout << " optimal sample = " << n << "\n";
    mk_statistic.reset();
    // starting value of variational parameters
    for (int i=0; i<varp.size(); ++i) 
      varp[i] = simulator.rng().random_real() * (varp_ub_[i]-varp_lb_[i]);
    //varp = simulator.varp_values();
    // Stochastic reconfiguration iterations
    for (unsigned iter=0; iter<sr_max_iter_; ++iter) {
      std::cout << " iter = " << iter << "\n";
      std::cout << " varp = " << varp.transpose() << "\n";
      double en = simulator.sr_function(varp, grad, sr_matrix);
      std::cout << " energy = " << en << "\n";
      std::cout << " grad = " << grad.transpose() << "\n";
      double gnorm = grad.squaredNorm();
      std::cout << " grad (sq) norm = " << gnorm << "\n";
      std::cout << " srmat = " << sr_matrix << "\n";
      // apply to stabiliser to sr matrix 
      for (int i=0; i<varp.size(); ++i) sr_matrix(i,i) += 1.0E-4;
      // search direction
      Eigen::VectorXd search_dir = sr_matrix.fullPivLu().solve(-grad);
      varp += sr_tstep_ * search_dir;
      // box constraint: project parameters into feasible region
      varp = varp_lb_.cwiseMax(varp.cwiseMin(varp_ub_));
      // add data to Mann-Kendall statistic
      mk_statistic << varp; 
      std::cout << " trend = " << mk_statistic.max_element_trend() << "\n"; 
      if (mk_statistic.is_full() && mk_statistic.max_element_trend()<sr_mktrend_tol_) {
        // converged
        mk_statistic.get_series_avg(varp);
        opt_varp_ << varp;
        break;
      }
    }
    // next sample
  }
  // print
  //opt_varp_.print_heading();
  //opt_varp_.print_result();

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

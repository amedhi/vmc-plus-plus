/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:19:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-09 17:46:04
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <string>
#include <stdexcept>
#include <boost/algorithm/string.hpp>
#include "./stochastic_reconf.h"

namespace vmc {

int StochasticReconf::init(const input::Parameters& inputs, const Simulator& simulator) 
{
  // problem size
  num_parms_ = simulator.num_varp();
  vparms_.resize(num_parms_);
  lbound_ = simulator.varp_lbound();
  ubound_ = simulator.varp_ubound();
  grad_.resize(num_parms_);
  sr_matrix_.resize(num_parms_,num_parms_);

  // optimization parameters
  int nowarn;
  num_sim_samples_ = inputs.set_value("measure_steps", 1000, nowarn);
  num_opt_samples_ = inputs.set_value("num_opt_samples", 30, nowarn);
  max_iter_ = inputs.set_value("sr_max_iter", 500, nowarn);
  start_tstep_ = inputs.set_value("sr_start_tstep", 0.05, nowarn);
  mk_series_len_ = inputs.set_value("sr_series_len", 40, nowarn);
  mk_thresold_ = inputs.set_value("sr_fluctuation_tol", 0.20, nowarn);
  refinement_cycle_ = max_iter_/5;
  if (mk_series_len_ > refinement_cycle_) {
    throw std::domain_error("StochasticReconf::init: sr_series_len > sr_max_iter/5");
  }

  // Mann-Kendall statistic
  mk_statistic_.resize(num_parms_);
  mk_statistic_.set_maxlen(mk_series_len_);

  // optimal parameter values
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  bool replace_mode = false;
  if (mode=="APPEND") replace_mode = false;
  optimal_parms_.init("OptParams", replace_mode);
  optimal_parms_.set_elements(simulator.varp_names());
  return 0;
}

int StochasticReconf::optimize(Simulator& simulator)
{
  // start optimization
  optimal_parms_.reset();
  for (unsigned n=0; n<num_opt_samples_; ++n) {
    std::cout << " optimal sample = " << n << "\n";
    // starting value of variational parameters
    vparms_ = (ubound_ - lbound_) * simulator.rng().random_real();
    // Stochastic reconfiguration iterations
    mk_statistic_.reset();
    for (unsigned iter=0; iter<max_iter_; ++iter) {
      std::cout << " iter = " << iter << "\n";
      std::cout << " varp = " << vparms_.transpose() << "\n";
      double en = simulator.sr_function(vparms_, grad_, sr_matrix_);
      std::cout << " energy = " << en << "\n";
      std::cout << " grad_ = " << grad_.transpose() << "\n";
      double gnorm = grad_.squaredNorm();
      std::cout << " grad_ (sq) norm = " << gnorm << "\n";
      std::cout << " srmat = " << sr_matrix_ << "\n";
      // apply to stabiliser to sr matrix 
      for (int i=0; i<num_parms_; ++i) sr_matrix_(i,i) += 1.0E-4;
      // search direction
      Eigen::VectorXd search_dir = sr_matrix_.fullPivLu().solve(-grad_);
      vparms_ += start_tstep_ * search_dir;
      // box constraint: project parameters into feasible region
      vparms_ = lbound_.cwiseMax(vparms_.cwiseMin(ubound_));
      // add data to Mann-Kendall statistic
      mk_statistic_ << vparms_; 
      std::cout << " trend = " << mk_statistic_.elem_max_trend() << "\n"; 
      if (mk_statistic_.is_full() && mk_statistic_.elem_max_trend()<mk_thresold_) {
        // converged
        mk_statistic_.get_series_avg(vparms_);
        optimal_parms_ << vparms_;
        break;
      }
    }
    // next sample
  }
  // print results
  //std::vector<std::string> xp({"x"});
  //optimal_parms_.print_heading(xp);
  std::vector<double> xv({simulator.hole_doping()});
  optimal_parms_.print_result(xv);
  return 0;
}



} // end namespace vmc
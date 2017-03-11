/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:19:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-10 22:26:27
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
  num_sim_samples_ = inputs.set_value("opt_measure_steps", 1000, nowarn);
  num_opt_samples_ = inputs.set_value("num_opt_samples", 30, nowarn);
  max_iter_ = inputs.set_value("sr_max_iter", 500, nowarn);
  start_tstep_ = inputs.set_value("sr_start_tstep", 0.05, nowarn);
  mk_series_len_ = inputs.set_value("sr_series_len", 40, nowarn);
  mk_thresold_ = inputs.set_value("sr_fluctuation_tol", 0.20, nowarn);
  grad_tol_ = inputs.set_value("sr_grad_tol", 0.05, nowarn);
  print_progress_ = inputs.set_value("sr_print_progress", true, nowarn);
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
  bool replace_mode = true;
  if (mode=="APPEND") replace_mode = false;
  optimal_parms_.init("OptParams", replace_mode);
  optimal_parms_.set_elements(simulator.varp_names());
  // observable file header
  std::stringstream heading;
  simulator.copyright_msg(heading);
  simulator.print_info(heading);
  std::vector<std::string> as_funct_of{"x"};
  optimal_parms_.print_heading(heading, as_funct_of);
  // progress file
  if (print_progress_) {
    progress_.open("log_optimization.txt");
    if (!progress_.is_open())
      throw std::runtime_error("StochasticReconf::init: file open failed");
    simulator.copyright_msg(progress_);
    simulator.print_info(progress_);
    progress_ << "#" << std::string(72, '-') << std::endl;
    progress_ << "Stochastic Reconfiguration" << std::endl;
    progress_ << "max_iter = " << max_iter_ << std::endl;
    progress_ << "start_tstep = " << start_tstep_ << std::endl;
    progress_ << "mk_series_len = " << mk_thresold_ << std::endl;
    progress_ << "grad_tol = " << grad_tol_ << std::endl;
    progress_ << "fluctuation_tol = " << mk_thresold_ << std::endl;
    progress_ << "optimization samples = " << num_opt_samples_ << std::endl;
    progress_ << "#" << std::string(72, '-') << std::endl;
  }
  return 0;
}

int StochasticReconf::optimize(Simulator& simulator)
{
  // start optimization
  optimal_parms_.reset();
  for (unsigned n=0; n<num_opt_samples_; ++n) {
    //std::cout << " optimal sample = " << n << "\n";
    if (print_progress_) {
      progress_ << "Starting sample " << n << " of " 
        << num_opt_samples_ << " ... " << std::flush;
    }
    // starting value of variational parameters
    vparms_ = lbound_+(ubound_-lbound_)*simulator.rng().random_real();
    // Stochastic reconfiguration iterations
    mk_statistic_.reset();
    double search_tstep = start_tstep_;
    int mc_samples = num_sim_samples_;
    unsigned iter;
    for (iter=0; iter<max_iter_; ++iter) {
      double en = simulator.sr_function(vparms_, grad_, sr_matrix_, mc_samples);
      // apply to stabiliser to sr matrix 
      for (unsigned i=0; i<num_parms_; ++i) sr_matrix_(i,i) += 1.0E-4;
      // search direction
      Eigen::VectorXd search_dir = sr_matrix_.fullPivLu().solve(-grad_);
      //getchar();
      // update variables
      vparms_ += search_tstep * search_dir;
      // box constraint and max_norm (of components not hitting boundary) 
      vparms_ = lbound_.cwiseMax(vparms_.cwiseMin(ubound_));
      /*double gnorm = std::abs(grad_[0]);
      for (unsigned i=0; i<num_parms_; ++i) {
        double x = vparms_(i);
        double lb = lbound_(i);
        double ub = ubound_(i);
        if (x < lb) vparms_(i) = lb;
        else if (x > ub) vparms_(i) = ub;
        else {
         double norm = std::abs(grad_[i]);
         if (max_norm < norm) max_norm = norm;
        }
      } */
      // add data to Mann-Kendall statistic
      //if (gnorm < grad_tol_) 
      mk_statistic_ << vparms_;
      double mk_trend = mk_statistic_.elem_max_trend();
      /*
      double gnorm = grad_.squaredNorm();
      std::cout << " iter = " << iter << "\n";
      std::cout << " search_dir = " << search_dir.transpose() << "\n";
      std::cout << " varp = " << vparms_.transpose() << "\n";
      std::cout << " energy = " << en << "\n";
      std::cout << " grad = " << grad_.transpose() << "\n";
      std::cout << " gnorm = " << gnorm << "\n";
      std::cout << " trend = " << mk_trend << "\n"; 
      */
      // convergence criteria
      if (mk_statistic_.is_full() && mk_trend<mk_thresold_) {
        // converged, add data point to store
        mk_statistic_.get_series_avg(vparms_);
        optimal_parms_ << vparms_;
        // print
        if (print_progress_) {
          double gnorm = grad_.squaredNorm();
          progress_ << "converged!" << std::endl;
          progress_ << "iteration = "<< iter << std::endl;
          progress_ << "MK trend = " << mk_statistic_.elem_max_trend() << "\n";
          progress_ << "gnorm = " << gnorm << std::endl;
          if (num_parms_<10) {
            progress_ << "varp = " << vparms_.transpose() << std::endl;
          }
          progress_ << "energy = " << en << std::endl;
          progress_ << std::endl;
        }
        break;
      }
      // refinement if not converged early
      if (iter % refinement_cycle_ == 0) {
        if (print_progress_) {
          progress_ << "next refinement cycle" << std::endl;
        }
        mc_samples *= 2;
        search_tstep *= 0.5;
      }
    }
    if (iter>=max_iter_ && print_progress_) {
      progress_ << "NOT converged):" << std::endl << std::endl;
    }
    // next sample
  }
  if (print_progress_) {
    if (optimal_parms_.num_samples()==num_opt_samples_)
      progress_ << "All done!" << std::endl;
    else if (optimal_parms_.num_samples()>0)
      progress_ << "Done. Some not converged." << std::endl;
    else
      progress_ << "None converged):" << std::endl;
  }
  progress_.close();
  // print results
  if (optimal_parms_.num_samples() > 0) {
    std::vector<double> xv({simulator.hole_doping()});
    optimal_parms_.print_result(xv);
    vparms_ = optimal_parms_.mean_data(); 
    return true;
  }
  else {
    std::cout << " StochasticReconf::optimize: optimization failed.\n";
    return false;
  }
}



} // end namespace vmc
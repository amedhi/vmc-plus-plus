/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-10 23:25:19
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "vmc.h"
#include <Eigen/SVD>

namespace vmc {

VMC::VMC(const input::Parameters& inputs) 
  : graph(inputs) 
  , model(inputs, graph.lattice()) 
  , config(inputs, graph, model)
  , num_varparms_{config.num_varparms()}
{
  // seed random generator
  config.rng().seed(inputs.set_value("rng_seed", 1));
  // mc parameters
  num_sites_ = graph.num_sites();
  num_measure_steps_ = inputs.set_value("measure_steps", 0);
  num_warmup_steps_ = inputs.set_value("warmup_steps", 0);
  min_interval_ = inputs.set_value("min_interval", 0);
  max_interval_ = inputs.set_value("max_interval", 0);
  check_interval_ = std::max(1,num_measure_steps_/10);

  // disorder
  if (site_disorder_.check(inputs)) {
    site_disorder_.init(inputs,graph,model,config,config.rng());
    model.add_disorder_term("disorder", model::op::ni_sigma());
    model.finalize(graph.lattice());
  }
  //if (model.have_disorder_term()) 

  // observables
  observables.init(inputs,copyright_msg,model,varp_names());
  observables.as_functions_of("x");

  if (observables.energy()) {
    term_energy_.resize(model.num_terms());
  }
  if (observables.energy_grad()) {
    grad_logpsi_.resize(num_varparms_);
    energy_grad_.resize(num_varparms_);
    energy_grad2_.resize(2*num_varparms_);
  } 

}

int VMC::disorder_start(const input::Parameters& inputs, 
  const unsigned& disorder_config, const bool& optimizing_mode, const bool& silent)
{
  site_disorder_.set_current_config(disorder_config);
  //site_disorder_.save_optimal_parms(config.vparm_values());
  bool with_psi_grad;
  if (observables.energy_grad()) with_psi_grad = true;
  else with_psi_grad = false;
  if (optimizing_mode) {
    observables.switch_off();
    observables.energy().switch_on();
  }
  silent_mode_ = silent;
  site_disorder_.get_optimal_parms();
  //return config.build(graph, inputs, with_psi_grad);
  return config.build(graph, site_disorder_.get_optimal_parms(), with_psi_grad);
}


int VMC::start(const input::Parameters& inputs, const bool& optimizing_mode, 
  const bool& silent)
{
  bool with_psi_grad;
  if (observables.energy_grad()) with_psi_grad = true;
  else with_psi_grad = false;
  if (optimizing_mode) {
    observables.switch_off();
    observables.energy().switch_on();
  }
  silent_mode_ = silent;
  return config.build(graph, inputs, with_psi_grad);
}

double VMC::energy_function(const Eigen::VectorXd& varp, Eigen::VectorXd& grad)
{
  // observables to calculate
  if (!observables.total_energy()) observables.total_energy().switch_on();
  bool with_psi_grad;
  if (grad.size()>0) {
    with_psi_grad = true;
    observables.energy_grad().switch_on(num_varparms_);
    observables.energy_grad2().switch_on(2*num_varparms_);
    grad_logpsi_.resize(num_varparms_);
  }
  else {
    observables.energy_grad().switch_off();
    with_psi_grad = false;
  }
  // build the config from the variational parameters
  config.build(graph, varp, with_psi_grad);
  run_simulation(observable_set::energy_grad);
  if (observables.energy_grad()) {
    finalize_energy_grad();
    grad = observables.energy_grad().mean_data();
  }
  //for (unsigned i=0; i<num_varparms_; ++i)
  //std::cout << " grad = " << grad.transpose() << "\n";
  //std::cout << " varp = " << x[0] << " " << x[1] << "\n";
  //std::cout << " energy = " << en << "\n";
  return observables.total_energy().mean();
}

double VMC::sr_function(const Eigen::VectorXd& varp, Eigen::VectorXd& grad, 
  Eigen::MatrixXd& sr_matrix, const int& sample_size)
{
  // observables to calculate
  if (!observables.total_energy()) observables.total_energy().switch_on();
  if (!observables.energy_grad()) {
    observables.energy_grad().switch_on(num_varparms_);
    observables.energy_grad2().switch_on(2*num_varparms_);
    grad_logpsi_.resize(num_varparms_);
    energy_grad_.resize(num_varparms_);
    energy_grad2_.resize(2*num_varparms_);
  } 
  if (!observables.sr_coeffs()) {
    observables.sr_coeffs().switch_on();
    // '\del(ln(psi))' plus upper triangular part of the sr_matrix 
    unsigned n = num_varparms_ + num_varparms_*(num_varparms_+1)/2;
    observables.sr_coeffs().set_elements(n);
    sr_coeffs_.resize(n);
  }
  // build the config from the variational parameters
  bool with_psi_grad = true;
  config.build(graph, varp, with_psi_grad);
  // run the simulation
  run_simulation(observable_set::sr_coeffs, sample_size);
  // gradient
  finalize_energy_grad();
  grad = observables.energy_grad().mean_data();
  // sr matrix
  sr_coeffs_ = observables.sr_coeffs().mean_data();
  unsigned k = num_varparms_;
  for (unsigned i=0; i<num_varparms_; ++i) {
    double x = sr_coeffs_[i];
    for (unsigned j=i; j<num_varparms_; ++j) {
      double y = sr_coeffs_[j];
      sr_matrix(i,j) = (sr_coeffs_[k] - x*y)/num_sites_;
      sr_matrix(j,i) = sr_matrix(i,j);
      ++k;
    }
  }
  return observables.total_energy().mean();
}

// simulation after optimization
int VMC::run_simulation(const Eigen::VectorXd& varp)
{
  observables.switch_off();
  observables.energy().switch_on();
  config.build(graph, varp);
  run_simulation();
  return 0;
}

int VMC::run_simulation(const observable_set& obs_set, const int& sample_size)
{
  // warmup run
  if (!silent_mode_) std::cout << " warming up... ";
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
  if (!silent_mode_) std::cout << "done\n";
  // measuring run
  int num_measurement = num_measure_steps_;
  if (sample_size>0) num_measurement = sample_size;
  int measurement_count = 0;
  int skip_count = min_interval_;
  observables.reset();
  while (measurement_count < num_measurement) {
    config.update_state();
    if (skip_count >= min_interval_) {
      if (config.accept_ratio()>0.5 || skip_count==max_interval_) {
        skip_count = 0;
        config.reset_accept_ratio();
        do_measurements(obs_set);
        ++measurement_count;
        if (!silent_mode_) print_progress(measurement_count);
      }
    }
    skip_count++;
  }
  if (!silent_mode_) {
    std::cout << " simulation done\n";
    config.print_stats();
  }
  return 0;
}


void VMC::print_results(void) 
{
  if (observables.energy_grad()) finalize_energy_grad();
  //observables.print(config.wavefunc().hole_doping());
  //observables.print_results(config.vparm_vector());
  observables.print_results(config.hole_doping());
  //std::cout << observables.energy().with_statistic() << "\n";
}

void VMC::print_progress(const int& num_measurement) const
{
  if (num_measurement%check_interval_==0)
  std::cout<<" measurement = "<< double(100.0*num_measurement)/num_measure_steps_<<" %\n";
}

/*
void VMC::warmup_config(void)
{
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
}
void VMC::update_config(void)
{
  int count = 0;
  config.reset_accept_ratio();
  while (count < max_interval_) {
    config.update_state();
    if (count>=min_interval_ && config.accept_ratio()>0.5) {
      config.reset_accept_ratio(); return;
    }
    ++count;
  }
}
*/

void VMC::init_observables(const input::Parameters& inputs)
{
  /*
  std::vector<std::string> elem_names;
  grad_logpsi_.resize(num_varparms_);
  energy_grad_.resize(num_varparms_);
  energy_grad2_.resize(2*num_varparms_);

  observables.init(inputs);
  model.get_term_names(elem_names);
  term_energy_.resize(model.num_terms());
  // energy 
  if (observables.energy()) {
    observables.energy().set_elements(elem_names);
  }
  // energy gradient
  if (observables.energy_grad()) {
    if (!observables.total_energy()) {
      observables.total_energy().switch_on();
      observables.total_energy().set_option(false);
    }
    observables.energy_grad2().switch_on();
    observables.energy_grad2().set_elements(2*num_varparms_);
    observables.energy_grad2().set_option(false);
    observables.energy_grad().set_elements(num_varparms_);
  }
  // in optimizing mode
  if (optimizing_mode_) {
    observables.switch_off();
    observables.energy().switch_on();
    observables.energy().set_elements(elem_names);
    observables.total_energy().switch_on();
    observables.energy_grad2().switch_on();
    observables.energy_grad2().set_elements(2*num_varparms_);
    observables.energy_grad2().set_option(false);
    observables.energy_grad().switch_on();
    observables.energy_grad().set_elements(num_varparms_);
  }

  need_energy_ = false;
  if (observables.energy()) need_energy_ = true;
  if (observables.energy_grad()) need_energy_ = true;
  if (observables.total_energy()) need_energy_ = true;
  if (observables.energy_grad()) need_gradient_ = true;
  else need_gradient_ = false;

  // open file & print heading
  observables.print_heading(config.vparm_names(), copyright_msg, model);
  */
}


void VMC::copyright_msg(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}



} // end namespace vmc






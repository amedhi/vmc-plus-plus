/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-08 21:11:09
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "simulator.h"
#include <Eigen/SVD>

namespace vmc {

Simulator::Simulator(const input::Parameters& inputs) 
  : graph(inputs) 
  , model(inputs, graph.lattice()) 
  , config(inputs, graph, model)
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
  // run mode
  optimizing_mode_ = inputs.set_value("optimizing_run",false);
  // observables
  observables.init(inputs,copyright_msg,model,varp_names());
  observables.as_functions_of("x");
  if (observables.energy()) term_energy_.resize(model.num_terms());
}

int Simulator::start(const input::Parameters& inputs, const bool& silent, 
  const bool& optimizing_mode)
{
  return config.build(inputs, graph, need_gradient_);
  silent_mode_ = silent;
  if (optimizing_mode) {
    observables.switch_off();
  }
}

double Simulator::energy_function(const var::parm_vector& x)
{
  Eigen::VectorXd v; 
  return energy_function(x, v);
}

double Simulator::energy_function(const var::parm_vector& x, Eigen::VectorXd& grad)
{
  bool need_gradient = false;
  //observables.energy_grad().switch_off();
  //if (grad.size()>0) {
    need_gradient = true;
    observables.energy_grad().switch_on();
 // }
  config.build(x, graph, need_gradient);
  run_simulation();
  if (need_gradient) {
    finalize_energy_grad();
    grad = observables.energy_grad().mean_data();
  }
  double en = observables.total_energy().mean();
  //for (unsigned i=0; i<num_varparms_; ++i)
  std::cout << " grad = " << grad.transpose() << "\n";
  std::cout << " varp = " << x[0] << " " << x[1] << "\n";
  std::cout << " energy = " << en << "\n";
  return en;
}

double Simulator::sr_function(const Eigen::VectorXd& vparms, Eigen::VectorXd& grad, 
  Eigen::MatrixXd& sr_matrix)
{
  observables.total_energy().switch_on();
  observables.energy_grad().switch_on();
  if (!observables.sr_matrix()) {
    observables.sr_matrix().switch_on();
    // upper triangular part of the sr_matrix plus '\del(ln(psi))' operators
    unsigned n = num_varparms_*(num_varparms_+1)/2+num_varparms_;
    observables.sr_matrix().set_elements(n);
    sr_matrix_el_.resize(n);
  }
  config.build(vparms, graph, true);
  run_simulation();
  finalize_energy_grad();
  // gradient
  grad.resize(num_varparms_);
  grad = observables.energy_grad().mean_data();
  // sr matrix
  sr_matrix.resize(num_varparms_, num_varparms_);
  sr_matrix_el_ = observables.sr_matrix().mean_data();
  unsigned k = num_varparms_;
  for (unsigned i=0; i<num_varparms_; ++i) {
    double x = sr_matrix_el_[i];
    for (unsigned j=i; j<num_varparms_; ++j) {
      double y = sr_matrix_el_[j];
      sr_matrix(i,j) = (sr_matrix_el_[k] - x*y)/num_sites_;
      sr_matrix(j,i) = sr_matrix(i,j);
      ++k;
    }
  }
  double en = observables.total_energy().mean();
  return en;
}

int Simulator::run_simulation(void)
{
  // run simulation
  int num_measurement = 0;
  int count = min_interval_;
  // warmup run
  if (!silent_mode_) std::cout << " warming up... ";
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
  if (!silent_mode_) std::cout << "done\n";
  // measuring run
  observables.reset();
  while (num_measurement < num_measure_steps_) {
    config.update_state();
    if (count >= min_interval_) {
      if (config.accept_ratio()>0.5 || count==max_interval_) {
        count = 0;
        config.reset_accept_ratio();
        do_measurements();
        ++num_measurement;
        if (!silent_mode_) print_progress(num_measurement);
      }
    }
    count++;
  }
  if (!silent_mode_) {
    std::cout << " simulation done\n";
    config.print_stats();
  }
  return 0;
}


void Simulator::print_results(void) 
{
  if (observables.energy_grad()) finalize_energy_grad();
  //observables.print(config.wavefunc().hole_doping());
  //observables.print_results(config.vparm_vector());
  observables.print_results(config.hole_doping());
  //std::cout << observables.energy().with_statistic() << "\n";
}

void Simulator::print_progress(const int& num_measurement) const
{
  if (num_measurement%check_interval_==0)
  std::cout<<" measurement = "<< double(100.0*num_measurement)/num_measure_steps_<<" %\n";
}

/*
void Simulator::warmup_config(void)
{
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
}
void Simulator::update_config(void)
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

void Simulator::init_observables(const input::Parameters& inputs)
{
  /*
  num_varparms_ = config.num_varparms();
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


void Simulator::copyright_msg(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}



} // end namespace vmc






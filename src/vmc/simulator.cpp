/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-25 13:41:51
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

  // no of variational parameters
  num_varparms_ = config.num_varparms();
  //std::cout << "num_varparms = " << num_varparms_ << "\n"; getchar();

  // observables
  observables.init(inputs);

  int nowarn;
  optimizing_mode_ = inputs.set_value("optimizing_run",false,nowarn);
  if (optimizing_mode_) {
    observables.switch_off();
    observables.energy().switch_on();
    grad_logpsi_.resize(num_varparms_);
    energy_grad_.resize(num_varparms_);
    energy_grad2_.resize(2*num_varparms_);
  }
  // sizes of vector observables
  std::vector<std::string> elem_names;
  model.get_term_names(elem_names);
  if (observables.energy()) {
    observables.energy().set_elements(elem_names);
  }
  if (observables.energy_grad()) {
    observables.total_energy().switch_on();
    num_varparms_ = config.num_varparms();
    grad_logpsi_.resize(num_varparms_);
    energy_grad_.resize(num_varparms_);
    energy_grad2_.resize(2*num_varparms_);
    observables.energy_grad().set_elements(2*num_varparms_);
  }
  need_energy_ = false;
  if (observables.energy()) need_energy_ = true;
  if (observables.energy_grad()) need_energy_ = true;
  if (observables.total_energy()) need_energy_ = true;
  if (need_energy_) config_energy_.resize(elem_names.size());

  // open file & print heading
  observables.print_heading(config.vparm_names(), copyright_msg, model);
}

int Simulator::init(const input::Parameters& inputs)
{
  return config.init(inputs, graph);
}

int Simulator::run(const bool& silent) 
{
  optimizing_mode_ = false;
  print_progress_ = !silent;
  if (observables.energy_grad()) {
    std::vector<double> varparms = config.vparm_values();
    int stat = config.init(varparms, graph, true);
    if (stat != 0) return -1;
  }
  run_simulation();
  if (observables.energy_grad()) finalize_energy_grad();
  return 0;
}

int Simulator::optimizing_run(const std::vector<double>& varparms, 
  const bool& need_energy_grad, const bool& silent)
{
  optimizing_mode_ = true;
  observables.total_energy().switch_on();
  if (need_energy_grad) {
    observables.energy_grad().switch_on();
    observables.energy_grad().set_elements(2*num_varparms_);
  }
  else observables.energy_grad().switch_off();
  int stat = config.init(varparms, graph, need_energy_grad);
  if (stat != 0) return -1;
  print_progress_ = !silent;
  run_simulation();
  if (need_energy_grad) finalize_energy_grad();
  return 0;
}

int Simulator::run_simulation(void)
{
  // run simulation
  int num_measurement = 0;
  int count = min_interval_;
  // warmup run
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
  if (print_progress_) std::cout << " warmup done\n";
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
        if (print_progress_) print_progress(num_measurement);
      }
    }
    count++;
  }
  if (print_progress_) {
    std::cout << " simulation done\n";
    config.print_stats();
  }
  return 0;
}

void Simulator::print_results(void) 
{
  //observables.print(config.wavefunc().hole_doping());
  observables.print(config.vparm_values());
  //std::cout << observables.energy().with_statistic() << "\n";
}

int Simulator::get_variational_parms(std::vector<double>& varparms)
{
  varparms = config.vparm_values();
  return varparms.size();
}

void Simulator::print_progress(const int& num_measurement) const
{
  if (num_measurement%check_interval_==0)
  std::cout<<" measurement = "<< double(100.0*num_measurement)/num_measure_steps_<<" %\n";
}

void Simulator::copyright_msg(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}



} // end namespace vmc






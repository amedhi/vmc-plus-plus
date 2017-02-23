/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-24 00:26:14
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "simulator.h"
#include <Eigen/SVD>

namespace vmc {

Simulator::Simulator(const input::Parameters& inputs) 
  : graph(inputs) 
  , model(inputs, graph.lattice()) 
  , config(inputs, graph, model)
  , observables_(inputs)
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

  // observables
  int nowarn;
  optimizing_mode_ = inputs.set_value("optimizing_run",false,nowarn);
  if (optimizing_mode_) {
    observables_.switch_off();
    observables_.energy().switch_on();
  }
  // sizes of vector observables
  if (observables_.energy()) {
    std::vector<std::string> elem_names;
    model.get_term_names(elem_names);
    config_energy_.resize(elem_names.size());
    observables_.energy().set_elements(elem_names);
  }
  observables_.print_heading(config.vparm_names(), copyright_msg, model);
}

int Simulator::init(const input::Parameters& inputs)
{
  return config.init(inputs, graph);
}

int Simulator::run(const bool& silent) 
{
  optimizing_mode_ = false;
  print_progress_ = !silent;
  run_simulation();
  return 0;
}

int Simulator::optimizing_run(const std::vector<double>& varparms, 
  const bool& need_energy_grad, const bool& silent)
{
  optimizing_mode_ = true;
  observables_.total_energy().switch_on();
  if (need_energy_grad) {
    observables_.energy_grad().switch_on();
    observables_.energy_grad().set_elements(config.vparm_names());
  }
  else observables_.energy_grad().switch_off();
  int stat = config.init(varparms, graph, need_energy_grad);
  if (stat != 0) return -1;
  print_progress_ = !silent;
  run_simulation();
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
  observables_.reset();
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
  //observables_.print(config.wavefunc().hole_doping());
  observables_.print(config.vparm_values());
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






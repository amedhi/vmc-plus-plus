/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-23 00:04:36
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "simulator.h"
#include <Eigen/SVD>

namespace vmc {

Simulator::Simulator(const input::Parameters& parms) 
  : graph(parms) 
  , model(parms, graph.lattice()) 
  , config(parms, graph, model)
{
  // seed random generator
  config.rng().seed(parms.set_value("rng_seed", 1));
  // mc parameters
  num_sites_ = graph.num_sites();
  num_measure_steps_ = parms.set_value("measure_steps", 0);
  num_warmup_steps_ = parms.set_value("warmup_steps", 0);
  min_interval_ = parms.set_value("min_interval", 0);
  max_interval_ = parms.set_value("max_interval", 0);
  check_interval_ = std::max(1,num_measure_steps_/10);

  // observables
  observables_.init(parms, model, print_copyright);
  //observables_.as_function_of("x");
  observables_.as_function_of(config.vparm_names());
  config_energy_.resize(model.num_terms());
}

int Simulator::run(const input::Parameters& inputs, const bool& silent) 
{
  int stat = config.init(inputs, graph);
  if (stat != 0) return -1;
  print_progress_ = !silent;
  run_simulation();
  return 0;
}

int Simulator::run(const std::vector<double>& varparms, const bool& silent) 
{
  int stat = config.init(varparms, graph);
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
  // output
  //observables_.print(config.wavefunc().hole_doping());
  observables_.print(config.vparm_values());
  return 0;
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


void Simulator::print_copyright(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}



} // end namespace vmc






/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-19 15:37:48
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

  // observables
  observables_.init(parms, model, print_copyright);
  observables_.as_function_of("x");
  config_energy_.resize(model.num_terms());
}

int Simulator::run(const input::Parameters& parms)
{
  int stat = config.init(parms);
  if (stat != 0) return -1;
  // run simulation
  int num_measurement = 0;
  int count = min_interval_;
  // warmup run
  for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
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
      }
    }
    count++;
  }
  // output
  observables_.print(0.0);
  std::cout << "total steps = " << config.num_updates() << "\n";
  return 0;
}


void Simulator::print_copyright(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}



} // end namespace vmc






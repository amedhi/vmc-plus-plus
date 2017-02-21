/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-21 22:56:42
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../mcdata/observables.h"
#include "./sysconfig.h"

namespace vmc {

class Simulator 
{
public:
  Simulator(const input::Parameters& parms); 
  ~Simulator() {}
  int init_config(void);
  int get_variational_parms(std::vector<double>& varparms);
  int run(const input::Parameters& parms);
  int run(const std::vector<double>& varparms);
  static void print_copyright(std::ostream& os);
private:
  lattice::LatticeGraph graph;
  model::Hamiltonian model;
  SysConfig config;
  unsigned num_sites_;

  // observables
  mc::Observables observables_;
  mutable mc::VectorData config_energy_;

  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_measure_steps_{0}; 
  int num_warmup_steps_{0};
  int min_interval_{0};
  int max_interval_{0};

  int run_simulation(void);
  int do_measurements(void);
  mc::VectorData config_energy(void) const;
};

} // end namespace vmc

#endif
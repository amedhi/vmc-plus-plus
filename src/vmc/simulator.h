/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-28 20:23:43
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
  Simulator(const input::Parameters& inputs); 
  ~Simulator() {}
  int start(const input::Parameters& inputs);
  int run(const bool& silent=false);
  int optimizing_run(const std::vector<double>& varparms, 
    const bool& need_energy_grad=false, const bool& silent=false);
  int get_variational_parms(std::vector<double>& varparms);
  const bool& optimizing_mode(void) const { return optimizing_mode_; }
  //void energy_gradient_off(void) { need_energy_grad_=false; }
  void print_results(void); 
  static void copyright_msg(std::ostream& os);
private:
  lattice::LatticeGraph graph;
  model::Hamiltonian model;
  SysConfig config;
  unsigned num_sites_;

  // optimization
  bool optimizing_mode_{false};

  // observables
  obs::ObservableSet observables;
  bool need_energy_{false};
  bool need_gradient_{false};
  mutable obs::vector config_energy_;
  mutable obs::vector energy_grad2_;
  mutable obs::vector energy_grad_;
  mutable RealVector grad_logpsi_;
  unsigned num_varparms_;

  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_measure_steps_{0}; 
  int num_warmup_steps_{0};
  int min_interval_{0};
  int max_interval_{0};
  int check_interval_{0};
  bool print_progress_{true};

  int run_simulation(void);
  int do_measurements(void);
  int finalize_energy_grad(void);
  void print_progress(const int& num_measurement) const;
  obs::vector config_energy(void) const;
};

} // end namespace vmc

#endif
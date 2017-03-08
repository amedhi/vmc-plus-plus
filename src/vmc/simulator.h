/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-08 21:11:04
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./observables.h"
#include "./sysconfig.h"

namespace vmc {

class Simulator 
{
public:
  Simulator(const input::Parameters& inputs); 
  ~Simulator() {}
  int start(const input::Parameters& inputs, const bool& silent=false, 
    const bool& optimizing_mode=false);
  int run_simulation(void);
  double energy_function(const var::parm_vector& x);
  double energy_function(const var::parm_vector& x, Eigen::VectorXd& grad);
  double sr_function(const Eigen::VectorXd& vparms, Eigen::VectorXd& grad, 
    Eigen::MatrixXd& sr_matrix);
  //void get_vparm_values(var::parm_vector& varparms) 
  //  { varparms = config.vparm_values(); }
  const unsigned& num_varparms(void) const { return config.num_varparms(); } 
  const var::parm_vector& varp_values(void) { return config.vparm_values(); }
  const var::parm_vector& varp_lbound(void) { return config.vparm_lbound(); }
  const var::parm_vector& varp_ubound(void) { return config.vparm_ubound(); }
  const std::vector<std::string>& varp_names(void) const { return config.vparm_names(); }
  RandomGenerator& rng(void) const { return config.rng(); }

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
  using vector_data = Observable::data_t;
  ObservableSet observables;
  bool need_energy_{false};
  bool need_gradient_{false};
  mutable vector_data term_energy_;
  mutable vector_data energy_grad2_;
  mutable vector_data energy_grad_;
  mutable RealVector grad_logpsi_;
  mutable vector_data sr_matrix_el_;
  unsigned num_varparms_;

  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_measure_steps_{0}; 
  int num_warmup_steps_{0};
  int min_interval_{0};
  int max_interval_{0};
  int check_interval_{0};
  bool silent_mode_{false};

  void init_observables(const input::Parameters& inputs);
  //void warmup_config(void);
  //void update_config(void);
  int do_measurements(void);
  int finalize_energy_grad(void);
  int finalize_sr_matrix(void);
  void print_progress(const int& num_measurement) const;
  vector_data config_energy(void) const;
};

} // end namespace vmc

#endif
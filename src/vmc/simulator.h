/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-09 17:08:24
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
  int start(const input::Parameters& inputs, const bool& optimizing_mode=false, 
    const bool& silent=false);
  int run_simulation(const observable_set& obs_set=observable_set::normal,  
    const int& sample_size=-1);
  double energy_function(const Eigen::VectorXd& varp, Eigen::VectorXd& grad);
  double sr_function(const Eigen::VectorXd& vparms, Eigen::VectorXd& grad, 
    Eigen::MatrixXd& sr_matrix);
  //void get_vparm_values(var::parm_vector& varparms) 
  //  { varparms = config.vparm_values(); }
  const unsigned& num_varp(void) const { return config.num_varparms(); } 
  const var::parm_vector& varp_values(void) { return config.vparm_values(); }
  const var::parm_vector& varp_lbound(void) const { return config.vparm_lbound(); }
  const var::parm_vector& varp_ubound(void) const { return config.vparm_ubound(); }
  const std::vector<std::string>& varp_names(void) const { return config.vparm_names(); }
  RandomGenerator& rng(void) const { return config.rng(); }
  const double& hole_doping(void) const { return config.hole_doping(); }
  void print_results(void); 
  static void copyright_msg(std::ostream& os);
private:
  lattice::LatticeGraph graph;
  model::Hamiltonian model;
  SysConfig config;
  unsigned num_sites_;
  unsigned num_varparms_;

  // observables
  enum class obs_set {normal, energy_grad, sr_matrix};
  obs_set observable_set_{obs_set::normal};
  using vector_data = Observable::data_t;
  ObservableSet observables;
  bool need_energy_{false};
  bool need_gradient_{false};
  mutable vector_data term_energy_;
  mutable vector_data energy_grad2_;
  mutable vector_data energy_grad_;
  mutable vector_data sr_coeffs_;
  mutable RealVector grad_logpsi_;

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
  void do_measurements(const observable_set& obs_set);
  int finalize_energy_grad(void);
  void measure_energy_grad(const double& total_en);
  int finalize_sr_matrix(void);
  void print_progress(const int& num_measurement) const;
  vector_data get_energy(void) const;
};

} // end namespace vmc

#endif
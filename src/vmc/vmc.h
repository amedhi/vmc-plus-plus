/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-10 12:40:05
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef VMC_H
#define VMC_H

#include "../scheduler/worker.h"
#include "../scheduler/mpi_comm.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./observables.h"
#include "./sysconfig.h"
#include "./disorder.h"

namespace vmc {

class VMC 
{
public:
  VMC(const input::Parameters& inputs); 
  ~VMC() {}
  int start(const input::Parameters& inputs, const bool& optimizing_mode=false, 
    const bool& silent=false);
  int run_simulation(const observable_set& obs_set=observable_set::normal,  
    const int& sample_size=-1);
  int run_simulation(const Eigen::VectorXd& varp);
  double energy_function(const Eigen::VectorXd& varp, Eigen::VectorXd& grad);
  double sr_function(const Eigen::VectorXd& vparms, Eigen::VectorXd& grad, 
    Eigen::MatrixXd& sr_matrix, const int& sample_size=-1);
  //void get_vparm_values(var::parm_vector& varparms) 
  //  { varparms = config.vparm_values(); }
  const unsigned& num_varp(void) const { return config.num_varparms(); } 
  const var::parm_vector& varp_values(void) { return config.vparm_values(); }
  const var::parm_vector& varp_lbound(void) const { return config.vparm_lbound(); }
  const var::parm_vector& varp_ubound(void) const { return config.vparm_ubound(); }
  const std::vector<std::string>& varp_names(void) const { return config.varp_names(); }
  RandomGenerator& rng(void) const { return config.rng(); }
  const double& hole_doping(void) const { return config.hole_doping(); }
  void print_results(void); 
  std::ostream& print_info(std::ostream& os) const { return model.print_info(os); }
  static void copyright_msg(std::ostream& os);
  const bool& disordered_system(void) const { return site_disorder_.exists(); }
  const unsigned& num_disorder_configs(void) const { return site_disorder_.num_configs(); }

  // disordered case
  int disorder_start(const input::Parameters& inputs, const unsigned& disorder_config, 
    const bool& optimizing_mode=false, const bool& silent=false);
  void save_optimal_parms(const var::parm_vector& optimal_parms) 
    { site_disorder_.save_optimal_parms(optimal_parms); }
  bool optimal_parms_exists(const unsigned& config) 
    { return site_disorder_.optimal_parms_exists(config); } 
  void set_disorder_config(const unsigned& config) 
    { site_disorder_.set_current_config(config); }
private:
  lattice::LatticeGraph graph;
  model::Hamiltonian model;
  SysConfig config;
  SiteDisorder site_disorder_;
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

  // pair correlation function
  using bond_iterator = lattice::LatticeGraph::bond_iterator;
  using bond_pair = std::pair<bond_iterator,bond_iterator>;
  std::vector<bond_pair> bond_pairs;


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
  void init_measurements(void);
  void do_measurements(const observable_set& obs_set);
  int finalize_energy_grad(void);
  void measure_energy_grad(const double& total_en);
  int finalize_sr_matrix(void);
  void print_progress(const int& num_measurement) const;
  vector_data get_energy(void) const;
};

} // end namespace vmc

#endif
/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-18 13:54:54
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-21 10:53:21
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SYSCONFIG_H
#define SYSCONFIG_H

#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../variational/wavefunction.h"
#include "../variational/projector.h"
#include "../variational/matrix.h"
#include "./basisstate.h"

namespace vmc {

constexpr double dratio_cutoff(void) { return 1.0E-8; } 
constexpr double gfactor_cutoff(void) { return 1.0E-8; } 

class SysConfig : public BasisState
{
public:
  SysConfig(const input::Parameters& parms, const lattice::graph::LatticeGraph& graph, 
    const model::Hamiltonian& model);
  ~SysConfig() {}
  int init(const input::Parameters& inputs, const lattice::graph::LatticeGraph& graph);
  int init(const std::vector<double>& vparms, const lattice::graph::LatticeGraph& graph);
  const std::vector<std::string>& vparm_names(void) const; 
  const std::vector<double>& vparm_values(void) const; 
  const std::vector<double>& vparm_lbounds(void) const; 
  const std::vector<double>& vparm_ubounds(void) const; 
  int update_state(void);
  double accept_ratio(void);
  void reset_accept_ratio(void);
  amplitude_t apply(const model::op::quantum_op& op, const unsigned& site_i, 
    const unsigned& site_j, const int& bc_phase) const;
  amplitude_t apply(const model::op::quantum_op& qn_op, const unsigned& site_i) const;
  int apply_niup_nidn(const unsigned& i) const;
  const int& num_updates(void) const { return num_updates_; }
  const var::Wavefunction& wavefunc(void) const { return wf; }
  //var::VariationalParms& var_parms(void) { return wf.var_parms(); }
private:
  var::Wavefunction wf;
  var::WavefunProjector projector;
  Matrix psi_mat;
  Matrix psi_inv;
  mutable ColVector psi_row;
  mutable RowVector psi_col;
  mutable RowVector inv_row;
  unsigned num_sites_;
  unsigned num_upspins_;
  unsigned num_dnspins_;

  // variational parameters
  mutable std::vector<std::string> vparm_names_;
  mutable std::vector<double> vparm_values_;
  mutable std::vector<double> vparm_lb_;
  mutable std::vector<double> vparm_ub_;
  unsigned num_projector_parms_{0};
  unsigned num_total_parms_{0};

  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_updates_{0};
  //int num_total_steps_{0};
  int num_uphop_moves_{0};
  int num_dnhop_moves_{0};
  int num_exchange_moves_{0};
  int refresh_cycle_{100};
  long num_proposed_moves_[move_t::end];
  long num_accepted_moves_[move_t::end];
  long last_proposed_moves_;
  long last_accepted_moves_;

  // helper methods
  int init_config(void);
  int set_run_parameters(void);
  int do_upspin_hop(void);
  int do_dnspin_hop(void);
  int do_spin_exchange(void);
  int inv_update_upspin(const int& upspin, const ColVector& psi_row, 
    const amplitude_t& det_ratio);
  int inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
    const amplitude_t& det_ratio);
  amplitude_t apply_upspin_hop(const unsigned& i, const unsigned& j,
    const int& bc_phase) const;
  amplitude_t apply_dnspin_hop(const unsigned& i, const unsigned& j,
    const int& bc_phase) const;
  amplitude_t apply_sisj_plus(const unsigned& i, const unsigned& j) const;
};

} // end namespace vmc

#endif
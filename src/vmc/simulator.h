/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-17 23:21:44
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../variational/wavefunction.h"
#include "../variational/projector.h"
#include "../variational/matrix.h"
#include "../mcdata/observables.h"
#include "./random.h"
#include "./basisstate.h"

namespace vmc {

constexpr double dratio_cutoff(void) { return 1.0E-8; } 

class Simulator 
{
public:
  Simulator(const input::Parameters& parms); 
  ~Simulator() {}
  int init_config(void);
  int run(const input::Parameters& parms);
  static void print_copyright(std::ostream& os);
private:
  lattice::graph::LatticeGraph graph;
  model::Hamiltonian model;
  var::Wavefunction wf;
  var::WavefunProjector projector;
  BasisState bstate;
  Matrix psi_mat;
  Matrix psi_inv;
  ColVector psi_row;
  RowVector psi_col;
  RowVector inv_row;
  unsigned num_sites_;
  unsigned num_upspins_;
  unsigned num_dnspins_;

  // observables
  mc::Observables observables_;
  //mc::VectorData energy_;

  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_measure_steps_{0}; 
  int num_warmup_steps_{0};
  int min_interval_{0};
  int max_interval_{0};
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
  int update_state(void);
  int set_run_parameters(void);
  double accept_ratio(void);
  void reset_accept_ratio(void);
  int do_upspin_hop(void);
  int do_dnspin_hop(void);
  int do_spin_exchange(void);
  int inv_update_upspin(const int& upspin, const ColVector& psi_row, 
    const amplitude_t& det_ratio);
  int inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
    const amplitude_t& det_ratio);
  int do_measurements(void);
};

} // end namespace vmc

#endif
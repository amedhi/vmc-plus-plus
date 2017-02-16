/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 12:52:07
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../variational/wavefunction.h"
#include "../variational/matrix.h"
#include "./random.h"
#include "./basisstate.h"

namespace vmc {

class Simulator 
{
public:
  Simulator(const input::Parameters& parms); 
  ~Simulator() {}
  int init_config(void);
  int run(const input::Parameters& parms);
private:
  lattice::graph::LatticeGraph graph;
  model::Hamiltonian model;
  var::Wavefunction wf;
  BasisState bstate;
  Matrix psi_mat;
  Matrix psi_inv;
  unsigned num_sites;
  unsigned num_upspins;
  unsigned num_dnspins;
  // mc parameters
  int num_data_samples; 
  int num_warmup_steps;
  int min_interval;

  // helper methods
  int update_state(void);
  int do_upspin_hop(void);
  int do_measurements(void);
};

} // end namespace simulator

#endif
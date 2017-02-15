/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef VMC_H
#define VMC_H

#include <iostream>
/*
#include <fstream>
#include <string>
#include <cmath>
#include <vector>
#include <array>
#include <Eigen/Core>
#include "../scheduler/worker.h"
//#include "./random.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../variational/wavefunction.h"
//#include "../observable/observables.h"
//#include "sitebasisstate.h"
//#include "observable_operator.h"
*/
#include "./simulator.h"

namespace vmc {

/*class Simulator : public scheduler::Worker
{
public:
  Simulator(input::Parameters& parms); 
  ~Simulator() {}
  int start(input::Parameters& parms) override;
  void run(void) override {} 
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
private:
  //var::Wavefunction psi_var_;
};
*/

class VMC : public scheduler::Worker
{
public:
  //using SiteState = SiteBasisState;
  //using state_idx = SiteBasisState::state_idx;
  //VMC() {}
  VMC(const input::Parameters& parms); 
  ~VMC() {}
  //using Model::update_parameters;
  int start(input::Parameters& parms) override;
  void run(void) override {} 
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);

private:
  //lattice::graph::LatticeGraph graph_;
  //model::Hamiltonian model_;
  Simulator simulator;


  //std::vector<SiteBasisState> state;
  
  // mc parameters
  int measure_samples; 
  int warmup;
  int min_interval;

  // observables
  /*Observables observables;
  mc::VectorData energy_terms;
  mc::VectorData energy_terms_sq;
  SiteObsOperator magn_op;
  SiteObsOperator strain_op;
  bool need_energy{false};
  bool need_magn{false};

  void init(void);
  void init_state_random(void);
  void init_boltzmann_table(void);
  void update_state(void) { update_state_metropolis(); }
  void update_parameters(input::Parameters& parms);
  void set_magn_op(const std::string& op, const std::string& site="i")
   { magn_op.init(basis(), op, site); }
  void set_strain_op(const std::string& op, const std::string& site="i")
   { strain_op.init(basis(), op, site); }
  virtual mc::VectorData get_energy(void) const;
  virtual double get_magnetization(void);
  virtual double get_potts_magnetization(void);
  virtual double get_strain(void);
  */
};


} // end namespace monte carlo

#endif

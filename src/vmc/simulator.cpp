/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 13:15:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "simulator.h"
#include <Eigen/SVD>

namespace vmc {

Simulator::Simulator(const input::Parameters& parms) 
  : graph(parms) 
  , model(parms, graph.lattice()) 
  , wf(parms, graph)
  , bstate(graph.num_sites(), model.double_occupancy())
{
  // seed random generator
  bstate.rng().seed(parms.set_value("rng_seed", 1));
  // mc parameters
  num_sites = graph.num_sites();
  num_data_samples = parms.set_value("data_samples", 0);
  num_warmup_steps = parms.set_value("warmup_steps", 0);
  min_interval = parms.set_value("min_interval", 0);
}

int Simulator::init_config(void) 
{
  num_upspins = wf.num_upspins();
  num_dnspins = wf.num_dnspins();
  if (num_sites==0) return -1;
  if (num_upspins==0 && num_dnspins==0) return -1;
  if (num_upspins != num_dnspins) 
    throw std::range_error("*Simulator::init: unequal UP & DN spin case not implemented");

  //bstate.init_spins(num_upspins, num_dnspins);
  bstate.init_spins(num_upspins, num_dnspins);

  // amplitude matrix & inv for the state
  psi_mat.resize(num_upspins, num_dnspins);
  psi_inv.resize(num_upspins, num_dnspins);
  // try for a well condictioned amplitude matrix
  double rcond = 0.0;
  int num_attempt = 0;
  while (rcond < 1.0E-30) {
    bstate.set_random();
    wf.get_amplitudes(psi_mat,bstate.upspin_sites(),bstate.dnspin_sites());
    // reciprocal conditioning number
    Eigen::JacobiSVD<Matrix> svd(psi_mat);
    // reciprocal cond. num = smallest eigenval/largest eigen val
    rcond = svd.singularValues()(svd.singularValues().size()-1)/svd.singularValues()(0);
    //std::cout << "rcondition number = "<< rcond << "\n";
    if (++num_attempt > 1000) {
      throw std::underflow_error("*Simulator::init: Ill conditioned wave function matrix.");
    }
  }
  //std::cout << psi_mat;
  //std::cout << bstate;
  // inverse amplitude matrix
  psi_inv = psi_mat.inverse();

  return 0;
}

int Simulator::run(const input::Parameters& parms)
{
  int stat = init_config();
  if (stat != 0) return -1;
  // run simulation
  int num_measurement = 0;
  int count = min_interval;
  // warmup run
  for (int n=0; n<num_warmup_steps; ++n) update_state();
  // measuring run
  while (num_measurement < num_data_samples) {
    update_state();
    if (count++ == min_interval) {
      count = 0;
      ++num_measurement;
      do_measurements();
    }
  }
  return 0;
}

int Simulator::update_state(void)
{
  for (int n=0; n<num_upspins; ++n) do_upspin_hop();
  //for (int n=0; n<num_sites; ++n) do_dnspin_hop();
  //for (int n=0; n<num_sites; ++n) do_spin_exchange();
  return 0;
}

int Simulator::do_upspin_hop(void)
{
  // int upspin = bstate.random_upspin_hop();
  int upspin, to_site;
  std::make_pair(upspin, to_site) = bstate.random_upspin_hop();
  if (to_site < 0) return 0; // valid move not found
  //unsigned p = rng.get_random_upspin();
  //unsigned s = rng.get_random_uphole();
  //bstate.random_upspin_hop();
  // if (x < bstate.rng().random_real()) {
  //}
  // else {
    //bstate.undo_last_move();
  //}

  return 0;
}

int Simulator::do_measurements(void)
{
  return 0;
}

} // end namespace simulator










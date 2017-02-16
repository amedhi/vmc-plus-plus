/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 23:56:29
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "simulator.h"
#include <Eigen/SVD>

namespace vmc {

Simulator::Simulator(const input::Parameters& parms) 
  : graph(parms) 
  , model(parms, graph.lattice()) 
  , wf(parms, graph)
  , projector(parms)
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
  psi_row.resize(num_dnspins);
  psi_col.resize(num_upspins);
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

int Simulator::update_state(void)
{
  for (int n=0; n<num_upspins; ++n) do_upspin_hop();
  for (int n=0; n<num_dnspins; ++n) do_dnspin_hop();
  //for (int n=0; n<num_sites; ++n) do_dnspin_hop();
  //for (int n=0; n<num_sites; ++n) do_spin_exchange();
  return 0;
}

int Simulator::do_upspin_hop(void)
{
  // int upspin = bstate.random_upspin_hop();
  int upspin, to_site;
  std::tie(upspin, to_site) = bstate.random_upspin_hop();
  //std::cout << "upspin, to_site = " << upspin << " " << to_site << "\n";
  if (to_site < 0) return 0; // valid move not found
  // new row for this move
  wf.get_amplitudes(psi_row, to_site, bstate.dnspin_sites());
  amplitude_t det_ratio = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  if (std::abs(det_ratio) < dratio_cutoff()) return 0; // for safety
  double proj_ratio = projector.gw_ratio(bstate.dblocc_increament());
  amplitude_t weight_ratio = det_ratio * proj_ratio;
  double transition_proby = std::norm(weight_ratio);
  //std::cout << "W = " << transition_proby << "\n";
  if (bstate.rng().random_real()<transition_proby) {
    bstate.accept_last_move();
    // update matrices
    psi_mat.row(upspin) = psi_row;
    amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
    for (int i=0; i<upspin; ++i) {
      amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
      psi_inv.col(i) -= beta * psi_inv.col(upspin);
    }
    for (int i=upspin+1; i<num_upspins; ++i) {
      amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
      psi_inv.col(i) -= beta * psi_inv.col(upspin);
    }
    psi_inv.col(upspin) *= ratio_inv;
  }
  //std::cout << "---upspin hop move done----\n";
  return 0;
}

int Simulator::do_dnspin_hop(void)
{
  // int upspin = bstate.random_upspin_hop();
  int dnspin, to_site;
  std::tie(dnspin, to_site) = bstate.random_dnspin_hop();
  if (to_site < 0) return 0; // valid move not found
  // new row for this move
  wf.get_amplitudes(psi_col, bstate.upspin_sites(), to_site);
  amplitude_t det_ratio = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
  if (std::abs(det_ratio) < dratio_cutoff()) return 0; // for safety
  double proj_ratio = projector.gw_ratio(bstate.dblocc_increament());
  amplitude_t weight_ratio = det_ratio * proj_ratio;
  double transition_proby = std::norm(weight_ratio);
  if (bstate.rng().random_real()<transition_proby) {
    bstate.accept_last_move();
    // update matrices
    psi_mat.col(dnspin) = psi_col;
    amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
    for (int i=0; i<dnspin; ++i) {
      amplitude_t beta = ratio_inv*psi_col.cwiseProduct(psi_inv.row(i)).sum();
      psi_inv.row(i) -= beta * psi_inv.row(dnspin);
    }
    for (int i=dnspin+1; i<num_upspins; ++i) {
      amplitude_t beta = ratio_inv*psi_col.cwiseProduct(psi_inv.row(i)).sum();
      psi_inv.row(i) -= beta * psi_inv.row(dnspin);
    }
    psi_inv.row(dnspin) *= ratio_inv;
  }
  return 0;
}

int Simulator::do_measurements(void)
{
  return 0;
}

} // end namespace simulator










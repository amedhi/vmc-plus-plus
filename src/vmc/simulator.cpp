/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-15 23:10:28
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "simulator.h"

namespace vmc {

Simulator::Simulator(const input::Parameters& parms) 
  : graph(parms) 
  , model(parms, graph.lattice()) 
  , wf(parms, graph)
  , basis_state(graph.num_sites(), model.double_occupancy())
{
  num_sites = graph.num_sites();
  rng.seed(parms.set_value("rng_seed", 1));
  rng.set_site_generator(0,num_sites-1);
  init_config();
}

int Simulator::init_config(void) 
{
  num_upspins = wf.num_upspins();
  num_dnspins = wf.num_dnspins();
  if (num_upspins==0 && num_dnspins==0) return 0;
  num_upholes = num_sites - num_upspins;
  num_dnholes = num_sites - num_dnspins;
  rng.set_upspin_generator(0,num_upspins-1);
  rng.set_dnspin_generator(0,num_dnspins-1);
  rng.set_uphole_generator(0,num_upholes-1);
  rng.set_dnhole_generator(0,num_dnholes-1);
  if (num_upspins != num_dnspins) {
    throw std::range_error("*Simulator::init: unequal UP & DN spin case not implemented");
  }
  // amplitude matrix & inv for the state
  psi_mat.resize(num_upspins, num_dnspins);
  psi_inv.resize(num_upspins, num_dnspins);

  basis_state.init_random(num_upspins, num_dnspins, rng);
  psi_mat = wf.amplitudes(basis_state.upspin_pos(),basis_state.dnspin_pos());
  std::cout << psi_mat;


  std::cout << basis_state;
  return 0;
}


} // end namespace simulator










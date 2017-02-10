/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-10 23:17:31
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "wavefunction.h"

namespace var {

Wavefunction::Wavefunction(const input::Parameters& inputs, 
  const lattice::graph::LatticeGraph& graph)
  : mf_model_(inputs, graph)
  , blochbasis_(graph)
  , num_kpoints_(blochbasis_.num_kpoints())
  , block_dim_(blochbasis_.subspace_dimension())
{
  if (mf_model_.is_pairing()) {
    bcs_init();
    if (block_dim_==1) type_ = wf_type::bcs_oneband;
    else type_ = wf_type::bcs_multiband;
  }
  else {
    type_ = wf_type::fermisea;
  }
  compute_amplitudes(graph);
}

//Wavefunction::update(const std::vector<double>& vparms)
//{
  // mf_model_.update(variational_parms);
  // compute_amlitudes()
//}

void Wavefunction::compute_amplitudes(const lattice::graph::LatticeGraph& graph)
{
  switch (type_) {
    case wf_type::bcs_oneband: 
      bcs_oneband(); 
      pair_amplitudes(graph);
      break;
    case wf_type::bcs_multiband: 
      bcs_multiband(); 
      pair_amplitudes(graph);
      break;
    case wf_type::fermisea: 
      fermisea(); 
      fermisea_amplitudes(graph); 
      break;
  }
  //(this->*construct_groundstate)();
}

void Wavefunction::pair_amplitudes(const lattice::graph::LatticeGraph& graph)
{
  unsigned num_sites = graph.num_sites();
  double one_by_nk = 1.0/static_cast<double>(num_kpoints_);
  for (unsigned i=0; i<num_sites; ++i) {
    unsigned m = graph.site_uid(i);
    auto Ri = graph.site_cellcord(i);
    for (unsigned j=0; j<num_sites; ++j) {
      unsigned n = graph.site_uid(j);
      auto Rj = graph.site_cellcord(j);
      std::complex<double> ksum(0.0);
      for (unsigned k=0; k<num_kpoints_; ++k) {
        Vector3d kvec = blochbasis_.kvector(k);
        ksum += cphi_k[k](m,n) * std::exp(ii()*kvec.dot(Ri-Rj));
      }
      psi_up_(i,j) = ksum * one_by_nk;
    }
  }
}


} // end namespace var












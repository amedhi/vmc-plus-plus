/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-17 21:40:13
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
  , num_sites_(graph.num_sites())
  , num_spins_(0)
  , num_upspins_(0)
  , num_dnspins_(0)
{
  set_particle_num(inputs);
  if (mf_model_.is_pairing()) {
    bcs_init(graph);
    if (block_dim_==1) type_ = wf_type::bcs_oneband;
    else type_ = wf_type::bcs_multiband;
  }
  else {
    type_ = wf_type::fermisea;
  }
  psi_up_.resize(num_sites_,num_sites_);
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
}

void Wavefunction::pair_amplitudes(const lattice::graph::LatticeGraph& graph)
{
  double one_by_nk = 1.0/static_cast<double>(num_kpoints_);
  for (unsigned i=0; i<num_sites_; ++i) {
    unsigned m = graph.site_uid(i);
    auto Ri = graph.site_cellcord(i);
    for (unsigned j=0; j<num_sites_; ++j) {
      unsigned n = graph.site_uid(j);
      auto Rj = graph.site_cellcord(j);
      std::complex<double> ksum(0.0);
      for (unsigned k=0; k<num_kpoints_; ++k) {
        Vector3d kvec = blochbasis_.kvector(k);
        ksum += cphi_k[k](m,n) * std::exp(ii()*kvec.dot(Ri-Rj));
      }
      psi_up_(i,j) = ksum * one_by_nk;
      //std::cout << psi_up_(i,j) << "\n"; 
      //getchar();
    }
  }
}

void Wavefunction::set_particle_num(const input::Parameters& inputs)
{
  band_filling_ = 1.0-inputs.set_value("x", 0.0);
  int num_sites = static_cast<int>(num_sites_);
  if (mf_model_.is_pairing()) {
    int n = static_cast<int>(std::round(0.5*band_filling_*num_sites));
    if (n<0 || n>num_sites) throw std::range_error("Wavefunction:: hole doping 'x' out-of-range");
    num_upspins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_upspins_;
    num_spins_ = num_upspins_ + num_dnspins_;
    band_filling_ = static_cast<double>(2*n)/num_sites;
  }
  else{
    int n = static_cast<int>(std::round(band_filling_*num_sites));
    if (n<0 || n>2*num_sites) throw std::range_error("Wavefunction:: hole doping 'x' out-of-range");
    num_spins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_spins_/2;
    num_upspins_ = num_spins_ - num_dnspins_;
  }
}

void Wavefunction::get_amplitudes(Matrix& psi, const std::vector<int>& row, 
  const std::vector<int>& col) const
{
  for (int i=0; i<row.size(); ++i)
    for (int j=0; j<col.size(); ++j)
      psi(i,j) = psi_up_(row[i],col[j]);
}

void Wavefunction::get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const
{
  for (int j=0; j<col.size(); ++j)
    psi_vec[j] = psi_up_(irow,col[j]);
}

void Wavefunction::get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const
{
  for (int j=0; j<row.size(); ++j)
    psi_vec[j] = psi_up_(row[j],icol);
}

void Wavefunction::get_amplitudes(amplitude_t& elem, const int& irow, 
  const int& jcol) const
{
  elem = psi_up_(irow,jcol);
}

} // end namespace var












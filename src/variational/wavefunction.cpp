/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-12 21:20:23
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "wavefunction.h"

namespace var {

Wavefunction::Wavefunction(const lattice::LatticeGraph& graph,
  const input::Parameters& inputs)
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
}

int Wavefunction::compute(const lattice::LatticeGraph& graph, 
  const input::Parameters& inputs, const bool& psi_gradient)
{
  set_particle_num(inputs);
  mf_model_.update(inputs,graph);
  if (mf_model_.need_noninteracting_mu()) {
    double mu = get_noninteracting_mu();
    //std::cout << "mu = " << mu << "\n"; getchar();
    mf_model_.update_mu(mu, graph); 
  }
  compute_amplitudes(psi_up_,graph);
  // psi gradients
  if (psi_gradient) {
    unsigned num_parms = mf_model_.varparms().size();
    var::parm_vector pvector(num_parms);
    //for (auto& p : mf_model_.varparms()) pvector.push_back(p.value());
    unsigned i = 0;
    for (auto& p : mf_model_.varparms()) {
     pvector[i] = p.value();
     ++i;
    }

    compute_gradients(graph, pvector);
    have_gradients_ = true;
  }
  else {
    have_gradients_ = false;
  }
  return 0;
}

int Wavefunction::compute(const lattice::LatticeGraph& graph, const var::parm_vector& pvector,
  const unsigned& start_pos, const bool& psi_gradient)
{
  mf_model_.update(pvector,start_pos,graph);
  compute_amplitudes(psi_up_,graph);
  // psi gradients
  if (psi_gradient) {
    compute_gradients(graph, pvector, start_pos);
    have_gradients_ = true;
  }
  else {
    have_gradients_ = false;
  }
  return 0;
}

int Wavefunction::compute_gradients(const lattice::LatticeGraph& graph, 
  const var::parm_vector& pvector, const unsigned& start_pos)
{
  // Gradient of amplitudes wrt the variational parameters 
  // by numerical differentiation (central defference formula)
  unsigned num_parm = mf_model_.varparms().size();
  psi_gradients_.resize(num_parm);
  work_mat.resize(num_sites_,num_sites_);
  //double scale = 0.005;
  unsigned i = 0;
  for (const auto& p : mf_model_.varparms()) {
    psi_gradients_[i].resize(num_sites_,num_sites_);
    //double h = scale * (p.ubound()-p.lbound());
    double h = p.diff_h();
    double inv_2h = 0.5/h;
    double x = pvector[start_pos+i];
    mf_model_.update(p.name(), x+h, graph);
    compute_amplitudes(psi_gradients_[i], graph);
    mf_model_.update(p.name(), x-h, graph);
    compute_amplitudes(work_mat, graph);
    // model to original state
    mf_model_.update(p.name(), x, graph);
    // derivative
    psi_gradients_[i] -= work_mat;
    psi_gradients_[i] *= inv_2h;
    ++i;
  }
  return 0;
}

int Wavefunction::compute_amplitudes(Matrix& psi_mat, const lattice::LatticeGraph& graph)
{
  switch (type_) {
    case wf_type::bcs_oneband: 
      bcs_oneband(); 
      pair_amplitudes(graph, psi_mat);
      break;
    case wf_type::bcs_multiband: 
      bcs_multiband(); 
      pair_amplitudes(graph, psi_mat);
      break;
    case wf_type::fermisea: 
      fermisea(); 
      fermisea_amplitudes(graph); 
      break;
  }
  return 0;
}

void Wavefunction::pair_amplitudes(const lattice::LatticeGraph& graph, Matrix& psi_mat)
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
      psi_mat(i,j) = ampl_part(ksum) * one_by_nk;
      //std::cout << psi_up_(i,j) << "\n"; 
      //getchar();
    }
  }
}

void Wavefunction::set_particle_num(const input::Parameters& inputs)
{
  hole_doping_ = inputs.set_value("hole_doping", 0.0);
  band_filling_ = 1.0-hole_doping_;
  int num_sites = static_cast<int>(num_sites_);
  if (mf_model_.is_pairing()) {
    int n = static_cast<int>(std::round(0.5*band_filling_*num_sites));
    if (n<0 || n>num_sites) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_upspins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_upspins_;
    num_spins_ = num_upspins_ + num_dnspins_;
    band_filling_ = static_cast<double>(2*n)/num_sites;
  }
  else{
    int n = static_cast<int>(std::round(band_filling_*num_sites));
    if (n<0 || n>2*num_sites) throw std::range_error("Wavefunction:: hole doping out-of-range");
    num_spins_ = static_cast<unsigned>(n);
    num_dnspins_ = num_spins_/2;
    num_upspins_ = num_spins_ - num_dnspins_;
    band_filling_ = static_cast<double>(n)/num_sites;
  }
  hole_doping_ = 1.0 - band_filling_;
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

void Wavefunction::get_gradients(Matrix& psi_grad, const int& n, 
  const std::vector<int>& row, const std::vector<int>& col) const
{
  if (!have_gradients_) 
    throw std::logic_error("Wavefunction::get_gradients: gradients were not computed");
  for (int i=0; i<row.size(); ++i)
    for (int j=0; j<col.size(); ++j)
      psi_grad(i,j) = psi_gradients_[n](row[i],col[j]);
}


void Wavefunction::get_vparm_names(std::vector<std::string>& vparm_names, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : mf_model_.varparms()) {
    vparm_names[start_pos+i] = p.name(); ++i;
  }
}

void Wavefunction::get_vparm_values(var::parm_vector& vparm_values, 
  unsigned start_pos)
{
  mf_model_.refresh_varparms(); 
  unsigned i = 0;
  for (auto& p : mf_model_.varparms()) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void Wavefunction::get_vparm_vector(std::vector<double>& vparm_values, 
  unsigned start_pos)
{
  mf_model_.refresh_varparms(); 
  unsigned i = 0;
  for (auto& p : mf_model_.varparms()) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void Wavefunction::get_vparm_lbound(var::parm_vector& vparm_lb, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : mf_model_.varparms()) {
    vparm_lb[start_pos+i] = p.lbound(); ++i;
  }
}

void Wavefunction::get_vparm_ubound(var::parm_vector& vparm_ub, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : mf_model_.varparms()) {
    vparm_ub[start_pos+i] = p.ubound(); ++i;
  }
}



} // end namespace var












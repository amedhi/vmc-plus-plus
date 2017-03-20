/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 23:06:41
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-21 01:20:35
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./bcs_state.h"

namespace var {

BCS_State::BCS_State(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph) 
  : GroundState(true)
{
  init(order_type, inputs, graph);
}

int BCS_State::init(const bcs& order_type, const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph)
{
  // sites
  num_sites_ = graph.num_sites();
  // particle number
  set_particle_num(inputs);
  // infinity limit
  large_number_ = 1.0E+2;

  // build MF Hamiltonian
  order_type_ = order_type;
  varparms_.clear();
  mf_model_.init(graph.lattice());
  std::string name;
  double defval, lb, ub;
  using namespace model;
  model::CouplingConstant cc;
  if (order_type_==bcs::swave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_bondterm(name="pairing", cc="delta_sc", op::pair_create());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else if (order_type_==bcs::dwave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_up());
    cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
    mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else {
    throw std::range_error("BCS_State::BCS_State: unidefined bcs order");
  }
  // chemical potential
  int info;
  mf_model_.add_parameter(name="mu", defval=0.0, inputs, info);
  if (info == 0) noninteracting_mu_ = false;
  else noninteracting_mu_ = true;
  if (inputs.set_value("mu_variational", false, info)) 
    varparms_.add("mu", defval=0.0, -2.0, +1.0);
  // finalize MF Hamiltonian
  mf_model_.finalize(graph);

  // bloch basis
  blochbasis_.construct(graph);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix();
  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  delta_k_.resize(kblock_dim_,kblock_dim_);
  dphi_k_.resize(kblock_dim_,kblock_dim_);
  phi_k_.resize(num_kpoints_);
  for (unsigned k=0; k<num_kpoints_; ++k) 
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
  return 0;
}

void BCS_State::get_wf_amplitudes(const input::Parameters& inputs, Matrix& psi) 
{
  update_parameters(inputs);
  // k-space pair amplitudes
  if (kblock_dim_==1) {
    pair_amplitudes_oneband(phi_k_);
  }
  else {
    pair_amplitudes_multiband(phi_k_);
  }
  // 'lattice space' pair amplitudes
  //pair_amplitudes_sitebasis(graph, phi_k, psi)

  /*double one_by_nk = 1.0/static_cast<double>(num_kpoints_);
  for (unsigned i=0; i<num_sites_; ++i) {
    unsigned m = graph.site_uid(i);
    auto Ri = graph.site_cellcord(i);
    for (unsigned j=0; j<num_sites_; ++j) {
      unsigned n = graph.site_uid(j);
      auto Rj = graph.site_cellcord(j);
      std::complex<double> ksum(0.0);
      for (unsigned k=0; k<num_kpoints_; ++k) {
        Vector3d kvec = blochbasis_.kvector(k);
        ksum += phi_k[k](m,n) * std::exp(ii()*kvec.dot(Ri-Rj));
      }
      psi(i,j) = ampl_part(ksum) * one_by_nk;
      //std::cout << psi_up_(i,j) << "\n"; 
      //getchar();
    }
  }*/
}

void BCS_State::pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
  // psi = FTU_ * PHI_K * conjugate(transpose(FTU_))
  // PHI_K is block diagonal (k-th block is phi_k) 
  for (unsigned i=0; i<num_kpoints_; ++i) {
    for (unsigned j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (unsigned k=0; k<num_kpoints_; ++k) {
        work_ += FTU_(i,k) * phi_k[j] * std::conj(FTU_(i,k));
      }
      // copy transformed block
      for (unsigned m=0; m<kblock_dim_; ++m) {
        for (unsigned n=0; n<kblock_dim_; ++n) 
          psi(i+m,j+n) = ampl_part(work_(m,n));
      }
    }
  }
}


void BCS_State::update_parameters(const input::Parameters& inputs)
{
  set_particle_num(inputs);
  mf_model_.update_parameters(inputs);
  if (noninteracting_mu_) {
    double mu_0 = get_noninteracting_mu();
    mf_model_.update_site_parameter("mu", mu_0);
  }
}

void BCS_State::pair_amplitudes_oneband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for one band system 
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    double ek = std::real(mf_model_.quadratic_spinup_block()(0,0)); 
    auto delta_k = mf_model_.pairing_part()(0,0);
    mf_model_.construct_kspace_block(-kvec);
    ek += std::real(mf_model_.quadratic_spinup_block()(0,0));;
    delta_k += mf_model_.pairing_part()(0,0);
    delta_k *= 0.5;
    double deltak_sq = std::norm(delta_k);
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (deltak_sq<1.0E-12 && ek<0.0) {
      phi_k[k](0,0) = large_number_ * std::exp(ii()*std::arg(delta_k));
    }
    else {
      phi_k[k](0,0) = 2.0*delta_k/ek_plus_Ek;
    }
    //std::cout << ek << "\n";
    //std::cout << delta_k << "\n";
    //std::cout << k << " " << cphi_k[k](0,0) << "\n";
    //getchar();
  }
}

void BCS_State::pair_amplitudes_multiband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for multi-band system (INTRABAND pairing only)
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_up.eigenvectors().conjugate();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    dphi_k_.setZero();
    for (unsigned i=0; i<kblock_dim_; ++i) {
      double ek = es_k_up.eigenvalues()[i] + es_minusk_up.eigenvalues()[i];
      double deltak_sq = std::norm(delta_k_(i,i));
      double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
      if (deltak_sq<1.0E-12 && ek<0.0) {
        dphi_k_(i,i) = large_number_ * std::exp(ii()*std::arg(delta_k_(i,i)));
      }
      else {
        dphi_k_(i,i) = 2.0*delta_k_(i,i)/ek_plus_Ek;
      }
    }
    // bcs ampitudes in original basis 
    for (unsigned i=0; i<kblock_dim_; ++i) 
      work_.col(i) = es_k_up.eigenvectors().col(i) * dphi_k_(i,i);
    phi_k[k] = work_ * es_minusk_up.eigenvectors().transpose();
    //std::cout << delta_k << "\n";
  } 
}



} // end namespace var
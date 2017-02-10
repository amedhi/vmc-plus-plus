/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-09 22:48:45
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-10 23:32:08
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "wavefunction.h"

namespace var {

void Wavefunction::bcs_init(void)
{
  // resizing
  mat_work.resize(block_dim_,block_dim_);
  mat_delta_k.resize(block_dim_,block_dim_);
  mat_dphi_k.resize(block_dim_,block_dim_);
  mat_dphi_k.resize(block_dim_,block_dim_);
  cphi_k.resize(num_kpoints_);
  for (unsigned k=0; k<num_kpoints_; ++k) cphi_k[k].resize(block_dim_,block_dim_);
  bcs_large_number_ = 1.0E+3;
}

void Wavefunction::bcs_oneband(void)
{
  // BCS wave function for one band system 
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    hk.compute(mf_model_.quadratic_spinup_block());
    double ek = hk.eigenvalues()[0];
    auto delta_k = mf_model_.pairing_part()(0,0);
    mf_model_.construct_kspace_block(-kvec);
    hminusk.compute(mf_model_.quadratic_spinup_block());
    ek += hminusk.eigenvalues()[0];
    delta_k += mf_model_.pairing_part()(0,0);
    double deltak_sq = std::norm(delta_k);
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (ek_plus_Ek>1.0E-12) {
      cphi_k[k](0,0) = delta_k/ek_plus_Ek;
    }
    else {
      cphi_k[k](0,0) = bcs_large_number_ * std::exp(ii()*std::arg(delta_k));
    }
  }
}

// Assumption: TR symmetry of the Hamiltonian
void Wavefunction::bcs_multiband(void)
{
  // BCS wave function for multiband system (INTRABAND pairing only)
  for (unsigned k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    hk.compute(mf_model_.quadratic_spinup_block());
    // pairing part 
    mat_delta_k = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    hminusk.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    mat_work = mat_delta_k + mf_model_.pairing_part().transpose();
    // transform pairing part
    mat_delta_k = hk.eigenvectors().adjoint() * mat_work * 
      hminusk.eigenvectors().conjugate();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    mat_dphi_k.setZero();
    for (unsigned i=0; i<block_dim_; ++i) {
      double ek = hk.eigenvalues()[i] + hminusk.eigenvalues()[i];
      double deltak_sq = std::norm(mat_delta_k(i,i));
      double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
      if (ek_plus_Ek>1.0E-12) {
        mat_dphi_k(i,i) = mat_delta_k(i,i)/ek_plus_Ek;
      }
      else {
        mat_dphi_k(i,i) = bcs_large_number_ * std::exp(ii()*std::arg(mat_delta_k(i,i)));
      }
    }
    // bcs ampitudes in original basis 
    cphi_k[k] = hk.eigenvectors() * mat_dphi_k.diagonal() * hminusk.eigenvectors().transpose();
    //std::cout << delta_k << "\n";
  } 
}

void Wavefunction::bcs_disordered(void)
{
  
}


} // end namespace var












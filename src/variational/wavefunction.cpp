/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-10 18:04:39
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
  , cphi_k(num_kpoints_)
  , delk_mat(block_dim_,block_dim_)
  , dphi_mat(block_dim_,block_dim_)
  , hk(block_dim_)
{
  //if (mf_model_.order()==mf_order::pairing) {
  //}
  // mf_model_.update(variational_parms);
  // for ()
  // mf_model_.construct_groundstate();
  for (unsigned k=0; k<num_kpoints_; ++k) cphi_k[k].resize(block_dim_,block_dim_);
  compute_amplitudes();
}

//Wavefunction::update(const std::vector<double>& vparms)
//{
  // mf_model_.update(variational_parms);
  // compute_amlitudes()
//}

void Wavefunction::compute_amplitudes(void)
{
  if (mf_model_.is_pairing()) {
    // BCS wave function
    for (unsigned k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      mf_model_.construct_kspace_block(kvec);
      hk.compute(mf_model_.quadratic_upspin_part());
      // transform pairing part
      delk_mat = hk.eigenvectors().adjoint() * mf_model_.pairing_part() * 
        hk.eigenvectors().conjugate();
      // 'dphi' amplitudes ignoring intraband pairing
      dphi_mat.setZero(); 
      for (unsigned i=0; i<block_dim_; ++i) {
        double ek = 2.0*hk.eigenvalues()[i];
        double deltak_sq = std::norm(delk_mat(i,i));
        double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
        dphi_mat(i,i) = 2.0*delk_mat(i,i)/ek_plus_Ek;
      }

      //mf_model_.kspace_transform(-kvec);
        //deltak_mat
        //phik_mat
      //std::cout << delta_k << "\n";
    } 
  }

  else {
  }
}

} // end namespace var












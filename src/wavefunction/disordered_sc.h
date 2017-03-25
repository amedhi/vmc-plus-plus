/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-22 22:41:54
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-25 15:25:29
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef DISORDERED_SC_H
#define DISORDERED_SC_H

#include <Eigen/Sparse>
#include "./groundstate.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class DisorderedSC : public GroundState
{
public:
  DisorderedSC() : GroundState(true) {}
  DisorderedSC(const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph); 
  ~DisorderedSC() {} 
  int init(const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph);
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_gradient(std::vector<Matrix>& psi_gradient) override; 
private:
  double large_number_{1.0E+2};
  unsigned num_bonds_;
  unsigned mu_start_;
  unsigned t_start_;
  unsigned delta_start_;
  // matrices
  using EigenTriplet = Eigen::Triplet<double>;
  std::vector<EigenTriplet> quadratic_coeffs_; 
  std::vector<EigenTriplet> pairing_coeffs_; 
  Eigen::SparseMatrix<double> quadratic_ham_;
  Eigen::SparseMatrix<double> pairing_ham_;
  Eigen::SparseMatrix<double> work_;
  RealMatrix delta_;
  Matrix psi_work_;
  std::vector<double> phi_k_;

  mutable Eigen::SelfAdjointEigenSolver<RealMatrix> es_quad_;
  mutable Eigen::SelfAdjointEigenSolver<RealMatrix> es_pair_;

  void get_pair_amplitudes_sitebasis(const Eigen::VectorXd& es_eigenvalues, 
  const RealMatrix& es_eigenvectors, const RealMatrix& delta, Matrix& psi);
};


} // end namespace var


#endif
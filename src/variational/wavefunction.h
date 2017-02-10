/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-10 17:40:04
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>
#include <map>
#include <stdexcept>
#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "mf_model.h"
#include "matrix.h"

namespace var {

// Wavefunction in 'Wannier space' (amplitudes in Wannier representation)
class Wavefunction 
{
public:
  //Wavefunction() {}
  Wavefunction(const input::Parameters& inputs, const lattice::graph::LatticeGraph& graph);
  ~Wavefunction() {}
  void update_amplitudes(const std::vector<double>& variatioanl_parm);
  void compute_grade(const std::vector<double>& variatioanl_parm);
private:
  MF_Model mf_model_;
  basis::BlochBasis blochbasis_;
  unsigned num_kpoints_;
  unsigned block_dim_;
  // BCS_state bcs_state_;
  // FS_state fermisea_;
  Matrix psi_up_;
  Matrix psi_dn_;

  // matrices & solvers
  std::vector<Matrix> cphi_k;
  Matrix delk_mat;
  Matrix dphi_mat;
  Matrix delta_k;
  Eigen::SelfAdjointEigenSolver<Matrix> hk;
  Eigen::SelfAdjointEigenSolver<Matrix> hminusk;
  Eigen::SelfAdjointEigenSolver<Matrix> hk_solver;
  Eigen::SelfAdjointEigenSolver<Matrix> bdg_solver;

  void compute_amplitudes(void);
};


} // end namespace var

#endif
/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-10 23:32:28
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

enum class wf_type {fermisea, bcs_oneband, bcs_multiband};

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
  wf_type type_;
  // BCS_state bcs_state_;
  // FS_state fermisea_;
  Matrix psi_up_;
  Matrix psi_dn_;

  // matrices & solvers
  double bcs_large_number_;
  Matrix mat_work;
  Matrix mat_delta_k;
  Matrix mat_dphi_k;
  std::vector<Matrix> cphi_k;
  Eigen::SelfAdjointEigenSolver<Matrix> hk;
  Eigen::SelfAdjointEigenSolver<Matrix> hminusk;

  void compute_amplitudes(const lattice::graph::LatticeGraph& graph);
  void pair_amplitudes(const lattice::graph::LatticeGraph& graph);
  void fermisea_amplitudes(const lattice::graph::LatticeGraph& graph) {}
  //void (Wavefunction::*construct_groundstate)(void);
  void bcs_init(void);
  void bcs_oneband(void);
  void bcs_multiband(void);
  void bcs_disordered(void);
  void fermisea(void) {}
};


} // end namespace var

#endif
/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-22 22:43:33
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
  Wavefunction(const lattice::LatticeGraph& graph, const input::Parameters& inputs);
  ~Wavefunction() {}
  void refresh_varparms(void) { mf_model_.refresh_varparms(); }
  const VariationalParms& varparms(void) const { return mf_model_.varparms(); }
  int compute(const lattice::LatticeGraph& graph, const input::Parameters& inputs,
    const bool& psi_gradient=false);
  int compute(const lattice::LatticeGraph& graph, const var::parm_vector& pvector,
    const unsigned& start_pos, const bool& psi_gradient=false);
  void compute_gradient(const std::vector<double>& vparms);
  const unsigned& num_upspins(void) const { return num_upspins_; }
  const unsigned& num_dnspins(void) const { return num_dnspins_; }
  const double& hole_doping(void) const { return hole_doping_; }
  void get_amplitudes(Matrix& psi, const std::vector<int>& row,  
    const std::vector<int>& col) const;
  void get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const;
  void get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const;
  void get_amplitudes(amplitude_t& elem, const int& irow, const int& jcol) const;
private:
  wf_type type_;
  MF_Model mf_model_;
  basis::BlochBasis blochbasis_;
  unsigned num_kpoints_;
  unsigned block_dim_;
  unsigned num_sites_;
  unsigned num_spins_;
  unsigned num_upspins_;
  unsigned num_dnspins_;
  double hole_doping_;
  double band_filling_;
  // BCS_state bcs_state_;
  // FS_state fermisea_;
  Matrix psi_up_;
  Matrix psi_dn_;
  Matrix work_mat;
  std::vector<Matrix> psi_gradients_;

  // matrices & solvers
  double bcs_large_number_;
  ComplexMatrix mat_work;
  ComplexMatrix mat_delta_k;
  ComplexMatrix mat_dphi_k;
  std::vector<ComplexMatrix> cphi_k;
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> hk;
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> hminusk;

  int compute_amplitudes(Matrix& psi_mat, const lattice::LatticeGraph& graph);
  void set_particle_num(const input::Parameters& inputs);
  double get_noninteracting_mu(void);
  void pair_amplitudes(const lattice::LatticeGraph& graph, Matrix& psi_mat);
  void fermisea_amplitudes(const lattice::LatticeGraph& graph) {}
  //void (Wavefunction::*construct_groundstate)(void);
  void bcs_init(const lattice::LatticeGraph& graph);
  void bcs_oneband(void);
  void bcs_multiband(void);
  void bcs_disordered(void);
  void fermisea(void) {}
};


} // end namespace var

#endif
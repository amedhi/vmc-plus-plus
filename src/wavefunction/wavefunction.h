/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-19 23:17:59
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>
//#include <map>
#include <memory>
#include <stdexcept>
#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "./mf_model.h"
#include "./matrix.h"
#include "./bcs_state.h"

namespace var {

enum class wf_type {normal, bcs_oneband, bcs_multiband, bcs_disordered};

// Wavefunction in 'Wannier space' (amplitudes in Wannier representation)
class Wavefunction 
{
public:
  //Wavefunction() {}
  Wavefunction(const lattice::LatticeGraph& graph, const input::Parameters& inputs,
    const bool& site_disorder=false);
  ~Wavefunction() {}
  const VariationalParms& varparms(void) const { return mf_model_.varparms(); }
  int compute(const lattice::LatticeGraph& graph, const input::Parameters& inputs, 
    const bool& psi_gradient=false);
  int compute(const lattice::LatticeGraph& graph, const var::parm_vector& pvector,
    const unsigned& start_pos, const bool& psi_gradient=false);
  //int compute_gradients(const lattice::LatticeGraph& graph);
  const unsigned& num_upspins(void) const { return num_upspins_; }
  const unsigned& num_dnspins(void) const { return num_dnspins_; }
  const double& hole_doping(void) const { return hole_doping_; }
  void get_vparm_names(std::vector<std::string>& names, unsigned start_pos) const; 
  void get_vparm_values(var::parm_vector& values, unsigned start_pos);
  void get_vparm_vector(std::vector<double>& vparm_values, unsigned start_pos);
  void get_vparm_lbound(var::parm_vector& lbounds, unsigned start_pos) const; 
  void get_vparm_ubound(var::parm_vector& ubounds, unsigned start_pos) const; 
  void get_amplitudes(Matrix& psi, const std::vector<int>& row,  
    const std::vector<int>& col) const;
  void get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const;
  void get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const;
  void get_amplitudes(amplitude_t& elem, const int& irow, const int& jcol) const;
  void get_gradients(Matrix& psi_grad, const int& n, 
    const std::vector<int>& row, const std::vector<int>& col) const;
private:
  std::unique_ptr<GroundState> ground_state_;
  //wf_descriptor wf_;
  std::string name_;
  bool pairing_type_{false};
  wf_type type_;
  basis::BlochBasis blochbasis_;
  MF_Model mf_model_;
  VariationalParms varparms_;
  bool need_noninteracting_mu_{true};
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
  bool have_gradients_{false};

  // matrices & solvers
  double bcs_large_number_;
  ComplexMatrix mat_work;
  ComplexMatrix mat_delta_k;
  ComplexMatrix mat_dphi_k;
  std::vector<ComplexMatrix> cphi_k;
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> hk;
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> hminusk;

  int compute_amplitudes(Matrix& psi_mat, const lattice::LatticeGraph& graph);
  int compute_gradients(const lattice::LatticeGraph& graph, 
    const var::parm_vector& pvector, const unsigned& start_pos=0);
  void set_particle_num(const input::Parameters& inputs);
  double get_noninteracting_mu(void);
  void pair_amplitudes(const lattice::LatticeGraph& graph, Matrix& psi_mat);
  void fermisea_amplitudes(const lattice::LatticeGraph& graph) {}
  //void (Wavefunction::*construct_groundstate)(void);
  void bcs_init(void);
  void bcs_oneband(void);
  void bcs_multiband(void);
  void bcs_disordered(const lattice::LatticeGraph& graph);
  void fermisea(void) {}
};


} // end namespace var

#endif
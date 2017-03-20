/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:32:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-20 23:31:20
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef GROUNDSTATE_H
#define GROUNDSTATE_H

#include <vector>
#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class GroundState
{
public:
  GroundState(const bool& pairing_type)
    : pairing_type_{pairing_type} {}
  ~GroundState() {} 
  virtual void get_wf_amplitudes(const input::Parameters& inputs, 
    Matrix& psi);
  const VariationalParms& varparms(void) { return varparms_; }
  const bool& pairing_type(void) const { return pairing_type_; }
  const unsigned& num_upspins(void) const { return num_upspins_; }
  const unsigned& num_dnspins(void) const { return num_dnspins_; }
  const double& hole_doping(void) const { return hole_doping_; }
protected:
  unsigned num_sites_{0};
  unsigned num_kpoints_{0};
  unsigned kblock_dim_{0};
  basis::BlochBasis blochbasis_;
  MF_Model mf_model_;
  VariationalParms varparms_;
  // fourier transform matrix
  ComplexMatrix FTU_;
  // solvers
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_minusk_up;

  void set_pairing_type(const bool& pairing_type) { pairing_type_=pairing_type; }
  void set_particle_num(const input::Parameters& inputs);
  void set_ft_matrix(void);
  double get_noninteracting_mu(void);
private:
  bool pairing_type_{false};
  unsigned num_spins_{0};
  unsigned num_upspins_{0};
  unsigned num_dnspins_{0};
  double hole_doping_{0.0};
  double band_filling_{1.0};
  //MF_Model mf_model_;
  //basis::BlochBasis blochbasis_;
};


} // end namespace var

#endif
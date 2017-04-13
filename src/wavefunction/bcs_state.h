/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:41:38
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-12 22:57:22
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef BCS_STATE_H
#define BCS_STATE_H

//#include <string>
//#include <complex>
//#include <vector>
//#include <map>
//#include <stdexcept>
//#include <Eigen/Eigenvalues>
//#include "../scheduler/task.h"
#include "./groundstate.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

enum class bcs {swave, dwave, af_swave, af_dwave};

class BCS_State : public GroundState
{
public:
  BCS_State() : GroundState(true) {}
  BCS_State(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph); 
  ~BCS_State() {} 
  int init(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph);
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_gradient(std::vector<Matrix>& psi_gradient) override; 
private:
  bcs order_type_;
  bool noninteracting_mu_{true};
  double large_number_{1.0E+2};
  // matrices
  ComplexMatrix work_;
  ComplexMatrix delta_k_;
  ComplexMatrix dphi_k_;
  std::vector<ComplexMatrix> phi_k_;
  std::vector<ComplexMatrix> work_k_;

  void get_pair_amplitudes_oneband(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_multiband(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, Matrix& psi);
  double get_mf_energy(void);
};


} // end namespace var


#endif
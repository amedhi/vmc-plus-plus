/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:41:38
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-20 01:01:54
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
  BCS_State() {}
  BCS_State(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph, const basis::BlochBasis& blochbasis); 
  ~BCS_State() {} 
  int init(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph, const basis::BlochBasis& blochbasis);
  int pair_amplitudes(std::vector<ComplexMatrix>& phi_k);
  bool pairing_type(void) const override { return true; }
private:
  bcs order_type_;
  MF_Model mf_model_;
  double large_number_{1.0E+2};
  // matrices
  ComplexMatrix work;
  ComplexMatrix delta_k;
  ComplexMatrix dphi_k;
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver_Hk;
  Eigen::SelfAdjointEigenSolver<ComplexMatrix> solver_Hmk;

};


} // end namespace var


#endif
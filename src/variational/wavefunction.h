/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   amedhi
* Last Modified time: 2017-01-30 22:26:55
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
#include "../scheduler/task.h"
#include "../model/model.h"
#include "matrix.h"

namespace var {

class Wavefunction 
{
public:
  Wavefunction() {}
  Wavefunction(const input::Parameters& inputs, const lattice::Lattice& lattice) 
  { construct(inputs, lattice); }
  ~Wavefunction() {}
  int construct(const input::Parameters& inputs, const lattice::Lattice& lattice);
private:
  using qn_op = model::qn_op;
  model::Model mf_model_;
  // 
  ComplexMatrix psi_up_;
  ComplexMatrix psi_dn_;
  int define_mf_model(const input::Parameters& inputs, const lattice::Lattice& lattice);
};


} // end namespace var

#endif
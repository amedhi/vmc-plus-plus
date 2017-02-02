/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-01 21:48:26
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
private:
  MF_Model mf_model_;
  Matrix psi_up_;
  Matrix psi_dn_;
};


} // end namespace var

#endif
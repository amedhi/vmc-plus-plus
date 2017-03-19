/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:32:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-20 00:55:20
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef GROUNDSTATE_H
#define GROUNDSTATE_H

#include <iostream>
//#include <string>
//#include <complex>
//#include <vector>
//#include <map>
//#include <stdexcept>
//#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class GroundState
{
public:
  //GroundState() {}
  virtual ~GroundState() {} 
  virtual bool pairing_type(void) const = 0;
  const VariationalParms& varparms(void) { return varparms_; }
  const unsigned& num_upspins(void) const { return num_upspins_; }
  const unsigned& num_dnspins(void) const { return num_dnspins_; }
protected:
  VariationalParms varparms_;
  unsigned num_sites_;
  unsigned num_spins_;
  unsigned num_upspins_;
  unsigned num_dnspins_;
  unsigned num_kpoints_;
  unsigned kblock_size_;
  //MF_Model mf_model_;
  //basis::BlochBasis blochbasis_;
};


} // end namespace var

#endif
/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-02 17:36:07
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MF_MODEL_H
#define MF_MODEL_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../model/model.h"
#include "../lattice/graph.h"
#include "blochbasis.h"
#include "matrix.h"

namespace var {
using name_value_pair = std::pair<std::string,double>;

class MF_Model : public model::Hamiltonian
{
public:
  MF_Model(const input::Parameters& inputs, const lattice::graph::LatticeGraph& graph);
  ~MF_Model() {}
  //void update_parameters(const var_parm& vparms_);
  void blochbasis_transform(const lattice::graph::LatticeGraph& graph);
private:
  using qn_op = model::qn_op;
  using Model = model::Hamiltonian;
  basis::BlochBasis kbasis_;
  std::vector<name_value_pair> vparms_;
  void make_variational(const std::vector<std::string>& pnames);
};


} // end namespace var

#endif
/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-09 22:20:27
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "wavefunction.h"

namespace var {

Wavefunction::Wavefunction(const input::Parameters& inputs, 
  const lattice::graph::LatticeGraph& graph)
  : mf_model_(inputs, graph)
{
  //if (mf_model_.order()==mf_order::pairing) {
  //}
  // mf_model_.update(variational_parms);
  // for ()
  // mf_model_.construct_groundstate();
  // compute_amlitudes()
}

//Wavefunction::update(const std::vector<double>& vparms)
//{
  // mf_model_.update(variational_parms);
  // compute_amlitudes()
//}
} // end namespace var












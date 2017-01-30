/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   amedhi
* Last Modified time: 2017-01-30 22:57:50
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "wavefunction.h"

namespace var {

int Wavefunction::construct(const input::Parameters& inputs, const lattice::Lattice& lattice) 
{
  define_mf_model(inputs, lattice);
  return 0;
}

int Wavefunction::define_mf_model(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  mf_model_.init(lattice);

  double defval;
  std::string name;
  model::CouplingConstant cc;
  
  // model parameters
  mf_model_.add_parameter(name="t", defval=1.0, inputs);
  mf_model_.add_parameter(name="delta_sc", defval=0.0, inputs);
  // bond operator terms
  mf_model_.add_bondterm(name="upspin_hop", cc="-t", qn_op::cdagicj_up);
  mf_model_.add_bondterm(name="dnspin_hop", cc="-t", qn_op::cdagicj_dn);
  mf_model_.add_bondterm(name="sc_pairing", cc="delta_sc", qn_op::cdagiup_cdagjdn);

  mf_model_.finalize(lattice);
  return 0;
}

} // end namespace var












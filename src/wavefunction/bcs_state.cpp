/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 23:06:41
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-20 01:02:36
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./bcs_state.h"

namespace var {

BCS_State::BCS_State(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph, const basis::BlochBasis& blochbasis) 
{
  init(order_type, inputs, graph, blochbasis);
}

int BCS_State::init(const bcs& order_type, const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph, const basis::BlochBasis& blochbasis)
{
  order_type_ = order_type;
  varparms_.clear();
  mf_model_.init(graph.lattice());

  std::string name;
  double defval, lb, ub;
  using namespace model;
  model::CouplingConstant cc;
  if (order_type_==bcs::swave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_bondterm(name="pairing", cc="delta_sc", op::pair_create());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else if (order_type_==bcs::dwave) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_up());
    cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
    mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    // variational parameters
    varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
  }
  else {
    throw std::range_error("BCS_State::BCS_State: unidefined bcs order");
  }
  // finalize MF Hamiltonian
  mf_model_.finalize(graph);

  // sizes
  num_kpoints_ = blochbasis.num_kpoints();
  kblock_size_ = blochbasis.subspace_dimension();
  work.resize(kblock_size_,kblock_size_);
  delta_k.resize(kblock_size_,kblock_size_);
  dphi_k.resize(kblock_size_,kblock_size_);
  large_number_ = 1.0E+2;

  return 0;
}

int BCS_State::pair_amplitudes(std::vector<ComplexMatrix>& phi_k)
{

  return 0;
}


} // end namespace var
/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-03 14:49:34
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "mf_model.h"

namespace var {

MF_Model::MF_Model(const input::Parameters& inputs, 
  const lattice::graph::LatticeGraph& graph)
  : kbasis_(graph)
{
  double defval;
  std::string name;
  model::CouplingConstant cc;
  
  Model::init(graph.lattice());
  // model parameters
  Model::add_parameter(name="mu", defval=0.0, inputs);
  Model::add_parameter(name="t", defval=1.0, inputs);
  Model::add_parameter(name="delta_sc", defval=0.0, inputs);
  // bond operator terms
  Model::add_bondterm(name="upspin_hop", cc="-t", qn_op::cdagicj_up);
  Model::add_bondterm(name="dnspin_hop", cc="-t", qn_op::cdagicj_dn);
  Model::add_bondterm(name="sc_pairing", cc="delta_sc", qn_op::cdagiup_cdagjdn);
  Model::finalize(graph.lattice());

  // variational parameters
  vparms_.clear();
  make_variational({"delta_sc", "mu"});
}

void MF_Model::make_variational(const std::vector<std::string>& pnames)
{
  for (const auto& pname : pnames) {
    vparms_.push_back({pname,get_parameter_value(pname)});
  }
}

void MF_Model::blochbasis_transform(const lattice::graph::LatticeGraph& graph)
{
  unsigned dim = kbasis_.block_dimension();
  Matrix ham(dim,dim);
  // bond terms
  unsigned i, j;
  lattice::graph::LatticeGraph::site_descriptor s, t;
  lattice::graph::LatticeGraph::out_edge_iterator ei, ei_end;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    switch (bterm->qn_operator()) {
      case qn_op::cdagicj_up:
      case qn_op::cdagicj_dn:
        for (i=0; i<dim; ++i) {
          s = kbasis_.state(i);
          for (std::tie(ei, ei_end)=graph.out_bonds(s); ei != ei_end; ++ei) {
            t = graph.target(ei);
            j = kbasis_.representative_state_idx(t);
            // mat(i,j) += bterm->coupling(graph.bond_type(ei));
          }
        }
      /*
        for (i=0; i<dim; ++i) {
          s = basis.site_basis(i);
          for (std::tie(ei, ei_end)=graph.out_edges(s); ei != ei_end; ++ei) {
            //std::cout << "edge type = " << graph.edge_type(ei) << "\n";
            t = graph.target_vertex(ei);
            j = basis.state_idx(t);
            if (j == basis.null_idx()) {
              type = graph.edge_type(ei);
              term = model.matrix_element(op, type);
              rs = basis.representative_state(t, graph, R);
              term = term * std::exp(+ii() * kvec.dot(R));
              j = basis.state_idx(rs);
              hk_up(i,j) += term;
              hk_up(j,i) += std::conj(term);
            }
          }
        }
      */
      default: break;
    }
  }

}

} // end namespace var












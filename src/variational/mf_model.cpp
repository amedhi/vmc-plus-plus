/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-18 07:01:13
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "mf_model.h"
#include <boost/algorithm/string.hpp>
//#include "../xml/pugixml.hpp"

namespace var {

MF_Model::MF_Model(const input::Parameters& inputs, 
  const lattice::graph::LatticeGraph& graph)
{
  // mean-field model
  define_model(inputs, graph);
  // 'unitcell representation' of the hamiltonian
  build_unitcell_terms(graph);

  // kspace matrices
  dim_ = graph.lattice().num_basis_sites();
  quadratic_block_up_.resize(dim_,dim_);
  pairing_block_.resize(dim_,dim_);
  work.resize(dim_,dim_);
  //work2.resize(dim_,dim_);
}

void MF_Model::define_model(const input::Parameters& inputs, const lattice::graph::LatticeGraph& graph)
{
  vparms_.clear();
  using namespace model;
  double defval;
  std::string name;
  model::CouplingConstant cc;
  // mean-field order & model
  Model::init(graph.lattice());
  std::string order_name = inputs.set_value("mf_order", "NONE");
  boost::to_upper(order_name);

  // chemical potential
  int info;
  add_parameter(name="mu", defval=0.0, inputs, info);
  if (info == 0) need_noninteracting_mu_ = false;
  else need_noninteracting_mu_ = true;
  if (inputs.set_value("mu_variational", false, info)) make_variational({"mu"});

  // assuming spin up-down symmetry, down-spin operators are not specified 
  if (order_name == "NONE") {
    order_ = mf_order::none;
    pairing_type_ = false;
    add_parameter(name="t", defval=1.0, inputs);
    add_bondterm(name="hopping", cc="-t", op::upspin_hop());
    add_siteterm(name="mu_term", cc="-mu", op::ni_up());
  }
  else if (order_name == "DWAVE_SC") {
    order_ = mf_order::dsc;
    pairing_type_ = true;
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="delta_sc", defval=1.0, inputs);
    add_bondterm(name="hopping", cc="-t", op::upspin_hop());
    add_siteterm(name="mu_term", cc="-mu", op::ni_up());
    cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
    add_bondterm(name="pairing", cc, op::pair_create());
    make_variational({"delta_sc"});
  }
  else if (order_name == "SWAVE_SC") {
    order_ = mf_order::ssc;
    pairing_type_ = true;
  }
  else if (order_name == "DISORDERED_SC") {
    order_ = mf_order::disorder_sc;
    pairing_type_ = true;
  }
  else {
    throw std::range_error("*error: mf_order: undefined order");
  }
  Model::finalize(graph.lattice());
}

void MF_Model::make_variational(const std::vector<std::string>& pnames)
{
  for (const auto& pname : pnames) {
    vparms_.push_back({pname,get_parameter_value(pname)});
  }
}

void MF_Model::build_unitcell_terms(const lattice::graph::LatticeGraph& graph)
{
  uc_siteterms_.resize(Model::num_siteterms());
  uc_bondterms_.resize(Model::num_bondterms());
  unsigned i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    uc_siteterms_[i].build_siteterm(*sterm, graph);
    i++;
  }
  i = 0;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    uc_bondterms_[i].build_bondterm(*bterm, graph);
    i++;
  }
}

void MF_Model::construct_kspace_block(const Vector3d& kvec)
{
  work.setZero(); 
  //work2 = Matrix::Zero(dim_,dim_);
  pairing_block_.setZero();
  // bond terms
  for (const auto& term : uc_bondterms_) {
    if (term.qn_operator().is_quadratic() && term.qn_operator().spin_up()) {
      for (unsigned i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        work += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
    if (term.qn_operator().is_pairing()) {
      for (unsigned i=0; i<term.num_out_bonds(); ++i) {
        Vector3d delta = term.bond_vector(i);
        pairing_block_ += term.coeff_matrix(i) * std::exp(ii()*kvec.dot(delta));
      }
    }
  }
  // add hermitian conjugate part
  quadratic_block_up_ = work + work.adjoint();
  // site terms 
  for (const auto& term : uc_siteterms_) {
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
    }
  }

  //quadratic_block_up_ += work1.adjoint();
  //pairing_block_ = work2;
  //pairing_block_ += work2.adjoint();
  // site terms
  //std::cout << "ek = " << quadratic_block_(0,0) << "\n";
}

void MF_Model::update_mu(const double& mu, const lattice::graph::LatticeGraph& graph)
{
  Model::update_parameter("mu", mu);
  unsigned i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    uc_siteterms_[i].build_siteterm(*sterm, graph);
    i++;
  }
}


/* Write a bond term like,
 H = \sum_{Ia,Jb}c^{\dag}_{Ia} t_{Ia,Jb} c_{Jb}
 for lattices with multiple sites per unit cell as
 H = \sum_{I,delta} Psi^{\dag}_{I} M^{\delta} Psi_{I+delta}
 Assumption: 'site's in the Graph are numbered contigously. For
 sites in the unitcell, the 'sl number' is same as 'uid'.
*/
void Unitcell_Term::build_bondterm(const model::HamiltonianTerm& hamterm,
  const lattice::graph::LatticeGraph& graph)
{
  unsigned dim = graph.lattice().num_basis_sites();
  lattice::graph::LatticeGraph::out_edge_iterator ei, ei_end;
  // get number of unique 'cell bond vectors'
  num_out_bonds_ = 0;
  for (unsigned i=0; i<dim; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      unsigned id = graph.vector_id(ei);
      if (id > num_out_bonds_) num_out_bonds_ = id;
    }
  }
  num_out_bonds_++;
  //std::cout << "num_out_bonds_ = " << num_out_bonds_ << "\n";
  bond_vectors_.resize(num_out_bonds_);
  coeff_matrices_.resize(num_out_bonds_);
  for (auto& M : coeff_matrices_) {
    M.resize(dim, dim);
    M = Matrix::Zero(dim,dim);
  }
  // operator
  op_ = hamterm.qn_operator();
  // build the matrices (for each 'bond vector')
  for (unsigned i=0; i<dim; ++i) {
    for (std::tie(ei, ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      unsigned id = graph.vector_id(ei);
      unsigned t = graph.target(ei);
      unsigned j = graph.site_uid(t);
      coeff_matrices_[id](i,j) += hamterm.coupling(graph.bond_type(ei));
      //std::cout << id << " " << coeff_matrices_[id](i,j) << "\n";
      bond_vectors_[id] = graph.vector(ei);
    }
  }
}

void Unitcell_Term::build_siteterm(const model::HamiltonianTerm& hamterm,
  const lattice::graph::LatticeGraph& graph)
{
  unsigned dim = graph.lattice().num_basis_sites();
  num_out_bonds_ = 0;
  bond_vectors_.resize(1);
  coeff_matrices_.resize(1);
  coeff_matrices_[0].resize(dim, dim);
  coeff_matrices_[0] = Matrix::Zero(dim,dim);
  // operator
  op_ = hamterm.qn_operator();
  // build the matrix 
  for (unsigned i=0; i<dim; ++i) {
    coeff_matrices_[0](i,i) = hamterm.coupling(graph.site_type(i));
  }
  bond_vectors_[0] = Vector3d(0,0,0);
}


/*void MF_Model::check_xml(void)
{
  std::cout << "Checking XML parser\n";
  pugi::xml_document doc;
  pugi::xml_parse_result result = doc.load_file("model.xml", pugi::parse_trim_pcdata);
  //std::cout << "Load result: " << result.description() << ", mesh name: " << "\n"; 
  pugi::xml_node model = doc.child("model");
  //std::cout << model.child("parameter").attribute("default").value() << std::endl;
  for (pugi::xml_node p = model.child("parameter"); p; p = p.next_sibling())
  {
    std::cout << "Parameter: ";
    for (pugi::xml_attribute attr = p.first_attribute(); attr; attr = attr.next_attribute())
    {
      std::cout << attr.name() << " = " << attr.as_double() << "\n";
    }
    std::cout << std::endl;
  }
  std::cout << "---------------------------------\n\n";
}*/



} // end namespace var












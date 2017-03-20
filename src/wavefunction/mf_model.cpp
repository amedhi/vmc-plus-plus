/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-20 16:36:23
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "mf_model.h"
#include <boost/algorithm/string.hpp>
#include "../expression/expression.h"
//#include "../xml/pugixml.hpp"

namespace var {

int MF_Model::finalize(const lattice::LatticeGraph& graph)
{
  Model::finalize(graph.lattice());
  dim_ = graph.lattice().num_basis_sites();
  quadratic_block_up_.resize(dim_,dim_);
  pairing_block_.resize(dim_,dim_);
  work.resize(dim_,dim_);
  build_unitcell_terms(graph);
  return 0;
}

MF_Model::MF_Model(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
{
  // mean-field model
  define_model(inputs, graph);
  // build 'unitcell representation' of the hamiltonian
  if (order_ == mf_order::disordered_sc) {
    dim_ = graph.num_sites();
  }
  else {
    dim_ = graph.lattice().num_basis_sites();
    build_unitcell_terms(graph);
  }
  // kspace block matrices
  // dim_ = blochbasis.subspace_dimension();
  quadratic_block_up_.resize(dim_,dim_);
  pairing_block_.resize(dim_,dim_);
  work.resize(dim_,dim_);
  //work2.resize(dim_,dim_);
}

void MF_Model::define_model(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
{
  using namespace model;
  double defval, lb, ub;
  std::string name;
  model::CouplingConstant cc;
  // mean-field order & model
  Model::init(graph.lattice());
  std::string order_name = inputs.set_value("mf_order", "NONE");
  boost::to_upper(order_name);

  // chemical potential term
  int info;
  add_parameter(name="mu", defval=0.0, inputs, info);
  if (info == 0) need_noninteracting_mu_ = false;
  else need_noninteracting_mu_ = true;
  if (inputs.set_value("mu_variational", false, info)) 
    make_variational("mu", -2.0, +1.0);

  // assuming spin up-down symmetry, down-spin operators are not specified 
  if (order_name == "NONE") {
    order_ = mf_order::none;
    pairing_type_ = false;
    add_parameter(name="t", defval=1.0, inputs);
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_siteterm(name="mu_term", cc="-mu", op::ni_up());
  }
  else if (order_name == "DWAVE_SC") {
    order_ = mf_order::dsc;
    pairing_type_ = true;
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="delta_sc", defval=1.0, inputs);
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_siteterm(name="mu_term", cc="-mu", op::ni_up());
    cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
    add_bondterm(name="pairing", cc, op::pair_create());
    // variational parameters
    make_variational("delta_sc", lb=0.0, ub=2.0);
  }
  else if (order_name == "SWAVE_SC") {
    order_ = mf_order::ssc;
    pairing_type_ = true;
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="delta_sc", defval=1.0, inputs);
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_bondterm(name="pairing", cc="delta_sc", op::pair_create());
    add_siteterm(name="mu_term", cc="-mu", op::ni_up());
    make_variational("delta_sc", lb=0.0, ub=2.0);
  }
  else if (order_name == "DISORDERED_SC") {
    order_ = mf_order::disordered_sc;
    pairing_type_ = true;
    add_bondterm(name="hopping", cc="-1.0", op::spin_hop());
    add_siteterm(name="site_mu", cc="0.0", op::ni_up());
    add_bondterm(name="pairing", cc="1.0", op::pair_create());
    // variational parameters
    need_noninteracting_mu_ = false;
    varparms_.clear();
    // hoping amplitudes
    bond_parms_start_ = 0;
    for (unsigned i=0; i<graph.num_bonds(); ++i) {
      name = "t_" + std::to_string(i);
      varparms_.add(name, -1.0, -2.0, 2.0);
    }
    // site potentials
    site_parms_start_ = varparms_.size();
    for (unsigned i=0; i<graph.num_sites(); ++i) {
      name = "mu_" + std::to_string(i);
      varparms_.add(name, 0.0, -4.0, 4.0);
    }
    // delta paraeters
    pair_parms_start_ = varparms_.size();
    for (unsigned i=0; i<graph.num_sites(); ++i) {
      name = "delta_" + std::to_string(i);
      varparms_.add(name, 0.5, -4.0, 4.0);
    }
  }
  else {
    throw std::range_error("*error: mf_order: undefined order");
  }
  Model::finalize(graph.lattice());
}

void MF_Model::update_parameters(const input::Parameters& inputs)
{
  Model::update_parameters(inputs);
  // update variational parameters
  for (auto& p : varparms_) 
    p.change_value(Model::get_parameter_value(p.name()));
  //build_unitcell_terms(graph);
  update_unitcell_terms();
}

void MF_Model::update(const parm_vector& pvector, const unsigned& start_pos, const lattice::LatticeGraph& graph)
{
  // for all variational parameters
  if (order_ == mf_order::disordered_sc) {
    unsigned i = 0;
    for (auto& p : varparms_) {
      p.change_value(pvector[start_pos+i]);
      ++i;
    }
    return;
  }

  unsigned i = 0;
  for (const auto& p : varparms_) {
    //std::cout << p.name() << " = " << pvector[start_pos+i] << "\n"; getchar();
    Model::update_parameter(p.name(), pvector[start_pos+i]);
    ++i;
  }
  //build_unitcell_terms(graph);
  update_unitcell_terms();
}

void MF_Model::update_parameter(const std::string& pname, const double& pvalue)
{
  Model::update_parameter(pname, pvalue, false);
  //std::cout << "updated ="<<pname<<"=" << get_parameter_value(pname) << "\n"; getchar();
  update_unitcell_terms();
}

void MF_Model::update_site_parameter(const std::string& pname, const double& pvalue)
{
  Model::update_parameter(pname, pvalue, false);
  for (unsigned i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
}

void MF_Model::make_variational(const std::string& name, const double& lb, const double& ub)
{
  varparms_.add(name, get_parameter_value(name), lb, ub);
}

void MF_Model::refresh_varparms(void) 
{
  for (auto& p : varparms_) p.change_value(get_parameter_value(p.name()));
}

void MF_Model::build_unitcell_terms(const lattice::LatticeGraph& graph)
{
  /*uc_siteterms_.resize(Model::num_siteterms());
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
  }*/
  usite_terms_.resize(Model::num_siteterms());
  ubond_terms_.resize(Model::num_bondterms());
  unsigned i = 0;
  for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
    usite_terms_[i].build_siteterm(*sterm, graph);
    i++;
  }
  i = 0;
  for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
    ubond_terms_[i].build_bondterm(*bterm, graph);
    i++;
  }
}

void MF_Model::construct_kspace_block(const Vector3d& kvec)
{
  work.setZero(); 
  pairing_block_.setZero();
  //work2 = Matrix::Zero(dim_,dim_);
  // bond terms
  //for (const auto& term : uc_bondterms_) {
  for (const auto& term : ubond_terms_) {
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
  //for (const auto& term : uc_siteterms_) {
  for (const auto& term : usite_terms_) {
    //std::cout << " --------- here --------\n";
    if (term.qn_operator().spin_up()) {
      quadratic_block_up_ += term.coeff_matrix();
      //std::cout << " sterm =" << term.coeff_matrix() << "\n"; //getchar();
    }
  }
  //quadratic_block_up_ += work1.adjoint();
  //pairing_block_ = work2;
  //pairing_block_ += work2.adjoint();
  // site terms
  //std::cout << "ek = " << quadratic_block_(0,0) << "\n";
}

const ComplexMatrix& MF_Model::quadratic_up_matrix(const lattice::LatticeGraph& graph) 
{
  if (num_bondterms() != 1)
    throw std::logic_error("MF_Model::quadratic_up_matrix: internal error");
  if (num_siteterms() != 1)
    throw std::logic_error("MF_Model::quadratic_up_matrix: internal error");

  work.setZero(); 
  int n = 0;
  // bond operator matrix elements
  for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int site_i = graph.source(b);
    int site_j = graph.target(b);
    // matrix elements 
    work(site_i, site_j) = varparms_[n].value();
    ++n;
  }
  quadratic_block_up_ = work + work.adjoint();
  // site operator matrix elements
  n = graph.num_bonds();
  for (unsigned i=0; i<graph.num_sites(); ++i) {
    quadratic_block_up_(i,i) += varparms_[n].value();
    ++n;
  }
  return quadratic_block_up_;
}

void MF_Model::get_pairing_varparms(ComplexMatrix& delta_k) const
{
  for (unsigned i=0; i<dim_; ++i)
    delta_k(i,i) = varparms_[pair_parms_start_+i].value();
}

void MF_Model::update_unitcell_terms(void)
{
  for (unsigned i=0; i<ubond_terms_.size(); ++i) 
    ubond_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
  for (unsigned i=0; i<usite_terms_.size(); ++i) 
    usite_terms_[i].eval_coupling_constant(Model::parameters(),Model::constants());
}


/* Write a bond term like,
 H = \sum_{Ia,Jb}c^{\dag}_{Ia} t_{Ia,Jb} c_{Jb}
 for lattices with multiple sites per unit cell as
 H = \sum_{I,delta} Psi^{\dag}_{I} M^{\delta} Psi_{I+delta}
 Assumption: 'site's in the Graph are numbered contigously. For
 sites in the unitcell, the 'sl number' is same as 'uid'.
*/

/*
void Unitcell_Term::build_bondterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  unsigned dim = graph.lattice().num_basis_sites();
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
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
    M.setZero();
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
  const lattice::LatticeGraph& graph)
{
  unsigned dim = graph.lattice().num_basis_sites();
  num_out_bonds_ = 0;
  bond_vectors_.resize(1);
  coeff_matrices_.resize(1);
  coeff_matrices_[0].resize(dim, dim);
  coeff_matrices_[0].setZero();
  // operator
  op_ = hamterm.qn_operator();
  // build the matrix 
  for (unsigned i=0; i<dim; ++i) {
    coeff_matrices_[0](i,i) = hamterm.coupling(graph.site_type(i));
  }
  bond_vectors_[0] = Vector3d(0,0,0);
}
*/

void UnitcellTerm::build_bondterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  dim_ = graph.lattice().num_basis_sites();
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  // get number of unique 'cell bond vectors'
  num_out_bonds_ = 0;
  for (unsigned i=0; i<dim_; ++i) {
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
    M.resize(dim_, dim_);
    M.setZero();
  }
  expr_matrices_.clear();
  expr_matrices_.resize(num_out_bonds_);
  for (auto& M : expr_matrices_) {
    M.resize(dim_);
    for (unsigned i=0; i<dim_; ++i) M[i].resize(dim_);
  }

  // operator
  op_ = hamterm.qn_operator();
  // build the matrices (for each 'bond vector')
  for (unsigned i=0; i<dim_; ++i) {
    for (std::tie(ei,ei_end)=graph.out_bonds(i); ei!=ei_end; ++ei) {
      unsigned id = graph.vector_id(ei);
      unsigned t = graph.target(ei);
      unsigned j = graph.site_uid(t);
      unsigned btype = graph.bond_type(ei);
      // expression
      std::string cc_expr(hamterm.coupling_expr(btype));
      boost::trim(cc_expr);
      if (cc_expr.size()>0) {
        if (cc_expr[0]!='-') expr_matrices_[id][i][j] += "+";
        expr_matrices_[id][i][j] += cc_expr;
      }
      // values
      coeff_matrices_[id](i,j) += hamterm.coupling(btype);
      //std::cout << id << " " << coeff_matrices_[id](i,j) << "\n";
      bond_vectors_[id] = graph.vector(ei);
    }
  }
}

void UnitcellTerm::build_siteterm(const model::HamiltonianTerm& hamterm,
  const lattice::LatticeGraph& graph)
{
  dim_ = graph.lattice().num_basis_sites();
  num_out_bonds_ = 1; // dummy, no real contrib
  bond_vectors_.resize(1);
  coeff_matrices_.resize(1);
  coeff_matrices_[0].resize(dim_,dim_);
  coeff_matrices_[0].setZero();
  expr_matrices_.resize(1);
  expr_matrices_[0].resize(dim_);
  for (unsigned i=0; i<dim_; ++i) expr_matrices_[0][i].resize(dim_);
  // operator
  op_ = hamterm.qn_operator();
  // build the matrix 
  for (unsigned i=0; i<dim_; ++i) {
    unsigned stype = graph.site_type(i);
    coeff_matrices_[0](i,i) = hamterm.coupling(stype);
    // expression
    std::string cc_expr(hamterm.coupling_expr(stype));
    boost::trim(cc_expr);
    if (cc_expr.size()>0) {
      if (cc_expr[0]!='-') expr_matrices_[0][i][i] += "+";
      expr_matrices_[0][i][i] += cc_expr;
    }
  }
  bond_vectors_[0] = Vector3d(0.0,0.0,0.0);
}

void UnitcellTerm::eval_coupling_constant(const model::ModelParams& pvals, const model::ModelParams& cvals)
{
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& p : pvals) {
    vars[p.first] = p.second;
    //std::cout << p.first << " = " << p.second << "\n"; getchar();
  }
  for (const auto& c : cvals) vars[c.first] = c.second;
  try { 
    for (unsigned n=0; n<num_out_bonds_; ++n) {
      for (unsigned i=0; i<dim_; ++i) {
        for (unsigned j=0; j<dim_; ++j) {
          std::string cc_expr(expr_matrices_[n][i][j]);
          if (cc_expr.size()>0) {
            coeff_matrices_[n](i,j) = expr.evaluate(cc_expr, vars); 
            //std::cout << "cc = " << coeff_matrices_[n](i,j) << "\n"; getchar();
          }
          else
            coeff_matrices_[n](i,j) = 0.0;
        }
      }
    }
  }
  catch (std::exception& e) 
  { 
    std::string msg = "UnitcellTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
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












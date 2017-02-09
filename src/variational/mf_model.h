/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-10 00:24:46
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
#include "../model/quantum_op.h"
#include "../model/model.h"
#include "../lattice/graph.h"
#include "blochbasis.h"
#include "matrix.h"

enum class mf_order {none, af, fm, ssc, dsc, pplusip, af_dsc, disorder_sc};

constexpr std::complex<double> ii(void) { return std::complex<double>{0.0,static_cast<double>(1.0)}; }

namespace var {
using name_value_pair = std::pair<std::string,double>;

//class SiteTerm : public std::unordered_map<unsigned, SiteOperatorTerm>
class Unitcell_Term
{
public:
  Unitcell_Term() {}
  ~Unitcell_Term() {}
  void build_bondterm(const model::HamiltonianTerm& sterm, const lattice::graph::LatticeGraph& graph);
  void build_siteterm(const model::HamiltonianTerm& sterm, const lattice::graph::LatticeGraph& graph);
  const unsigned& num_out_bonds(void) const { return num_out_bonds_; } 
  const Vector3d& bond_vector(const unsigned& i) const { return bond_vectors_[i]; }
  const Matrix& coeff_matrix(const unsigned& i=0) const { return coeff_matrices_[i]; }
  //const double& coupling(const unsigned& site_type) const; 
  const model::op::quantum_op& qn_operator(void) const { return op_; }
private:
  model::op::quantum_op op_;
  unsigned num_out_bonds_;
  std::vector<Matrix> coeff_matrices_;
  std::vector<Vector3d> bond_vectors_;
};

class MF_Model : public model::Hamiltonian
{
public:
  MF_Model(const input::Parameters& inputs, const lattice::graph::LatticeGraph& graph);
  ~MF_Model() {}
  //void update_parameters(const var_parm& vparms_);
  const bool& is_pairing(void) const { return pairing_type_; }
  void kspace_transorm(const unsigned& k);
  const Matrix& quadratic_block(void) const { return quadratic_block_; }
  const Matrix& pairing_block(void) const { return pairing_block_; }
private:
  using Model = model::Hamiltonian;
  mf_order order_;
  bool pairing_type_;
  basis::BlochBasis blochbasis_;
  std::vector<name_value_pair> vparms_;
  std::vector<Unitcell_Term> uc_siteterms_;
  std::vector<Unitcell_Term> uc_bondterms_;
  // matrices in kspace representation
  unsigned dim_;
  Matrix quadratic_block_;
  Matrix pairing_block_;
  Matrix work1, work2;

  //void check_xml(void);
  void define_model(const input::Parameters& inputs, const lattice::graph::LatticeGraph& graph);
  void deine_pairing(const std::vector<std::string>& pnames);
  void make_variational(const std::vector<std::string>& pnames);
  void build_unitcell_terms(const lattice::graph::LatticeGraph& graph);
};


} // end namespace var

#endif
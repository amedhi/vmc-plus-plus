/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-05 13:24:26
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
  const Matrix& coeff_matrix(const unsigned& i) const { return coeff_matrices_[i]; }
  //const double& coupling(const unsigned& site_type) const; 
  //const std::string& name(void) const { return name_; }
private:
  model::qn_op op_;
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
  void build_blochbasis_groundstate(void);
  //void blochbasis_transform(const lattice::graph::LatticeGraph& graph);
private:
  using qn_op = model::qn_op;
  using Model = model::Hamiltonian;
  std::vector<name_value_pair> vparms_;
  std::vector<Unitcell_Term> uc_siteterms_;
  std::vector<Unitcell_Term> uc_bondterms_;
  basis::BlochBasis blochbasis_;

  void make_variational(const std::vector<std::string>& pnames);
  void build_unitcell_terms(const lattice::graph::LatticeGraph& graph);
};


} // end namespace var

#endif
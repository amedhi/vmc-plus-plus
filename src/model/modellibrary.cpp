/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-13 00:18:03
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Hamiltonian::define_model(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  //int info;
  //unsigned ntypes;
  //std::vector<MatrixElement> matrix_elem(20);
  double defval;
  //unsigned sitetype, change, type, src_type, tgt_type;
  std::string name; //, matrixelem, op, qn, site, src, tgt, fact;
  //SiteBasis site_basis;
  //BasisDescriptor basis;
  //QuantumNumber::value_type min, max, step;
  CouplingConstant cc;

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    // model parameters
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="U", defval=0.0, inputs);
    // bond operator terms
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_siteterm(name="hubbard", cc="U", op::hubbard_int());
  }

  else if (model_name == "TJ") {
    mid = model_id::TJ;
    set_no_dbloccupancy();
    // model parameters
    add_parameter(name="t", defval=1.0, inputs);
    add_parameter(name="J", defval=0.0, inputs);
    // bond operator terms
    add_bondterm(name="hopping", cc="-t", op::spin_hop());
    add_bondterm(name="exchange", cc="J", op::sisj_plus());
  }

  /*------------- undefined lattice--------------*/
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  // if the model has disorder
  int nowarn;
  if (inputs.set_value("have_disorder",false,nowarn)) {
    set_have_disorder();
    add_siteterm(name="disorder", cc="0.0", op::ni_sigma());
    // file containing disorder potential
    // model signature string (to be used a folder name)
    std::ostringstream signature;
    signature << "L_" << lattice.size1() << "x" 
      << lattice.size2() << "x" << lattice.size3();
    signature << "_" << model_name;
    signature.precision(3);
    signature.setf(std::ios_base::fixed);
    for (const auto& p : parms_) info_str_ << "_" << p.first << p.second;
    disorder_file_ = inputs.set_value("input_prefix","system"); 
    disorder_file_ += "/";
    disorder_file_ += signature.str();
    disorder_file_ += "/";
    disorder_file_ += "disorder_potential.txt";
  }
  return 0;
}

int Hamiltonian::construct(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  init(lattice);
  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}

int Hamiltonian::init(const lattice::Lattice& lattice)
{
  // reset
  parms_.clear();
  //operators.clear();
  // maps of site & bond type values (to contigous type values)
  sitetypes_map_ = lattice.sitetypes_map();
  bondtypes_map_ = lattice.bondtypes_map();
  // maps of a given bond type to the types of its target
  /*bond_sites_map_.clear();
  for (unsigned i=0; i<lattice.num_basis_bonds(); ++i) {
    lattice::Bond b = lattice.basis_bond(i);
    lattice::Site src = lattice.basis_site(b.src_id());
    lattice::Site tgt = lattice.basis_site(b.tgt_id());
    bond_sites_map_.insert({b.type(), std::make_pair(src.type(), tgt.type())});
    //std::cout << "bond_site_map = "<<b.type()<<" "<<src.type()<<" "<<tgt.type()<<"\n";
  }*/
  return 0;
}


} // end namespace model

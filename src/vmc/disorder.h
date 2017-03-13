/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-13 11:22:16
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-13 22:41:59
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef DISORDER_H
#define DISORDER_H

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <vector>
#include "../lattice/graph.h"
#include "../model/model.h"

namespace vmc {

class SiteDisorder
{
public:
  SiteDisorder() {}
  SiteDisorder(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model) { init(inputs, graph, model); }
  ~SiteDisorder() {}
  int init(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model); 
  const unsigned& num_configs(void) const { return num_configs_; }
private:
  unsigned num_sites_;
  unsigned num_configs_;
  double bandwidth_{1.0};
  bool potential_set_{false};
  std::vector<std::vector<double> > disorder_pot_;
  std::string datafile_prefix_;
};





} // end namespace vmc

#endif
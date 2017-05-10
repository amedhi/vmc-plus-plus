/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:41:40
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-10 23:42:14
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_ENERGY_H
#define OBS_ENERGY_H

#include "../mcdata/mc_observable.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./sysconfig.h"

namespace vmc {

class Energy : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const model::Hamiltonian& model);
  void measure(const lattice::LatticeGraph& graph, const model::Hamiltonian& model,
    const SysConfig& config);
  const mcdata::data_t& config_value(void) const { return config_value_; }
private:
  mcdata::data_t config_value_;
};

class EnergyGradient : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const SysConfig& config);
  void measure(const SysConfig& config, const double& config_energy);
  void finalize(const double& mean_energy);
private:
  unsigned num_varp_{0};
  mcdata::data_t config_value_;
  RealVector grad_logpsi_;
  mcdata::MC_Observable grad_terms_;
};


} // end namespace vmc

#endif

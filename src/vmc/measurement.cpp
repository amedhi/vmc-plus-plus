/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-17 23:30:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-24 00:29:39
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include "simulator.h"

namespace vmc {

int Simulator::do_measurements(void)
{
  if (optimizing_mode_) {
    double energy = config_energy().sum();
    observables_.total_energy() << energy;
    if (observables_.energy_grad()) {
      //config.grad_log_psi()
    }
    return 0;
  }

  // normal run
  if (observables_.energy()) {
    observables_.energy() << config_energy();
  }
  return 0;
}

mc::VectorData Simulator::config_energy(void) const
{
  using op_id = model::op_id;
  //for (auto& elem : config_energy_) elem = 0.0;
  config_energy_.setZero();
  // bond energies
  if (model.has_bondterm()) {
    Matrix matrix_elem(model.num_bondterms(),graph.num_bond_types());
    matrix_elem.setZero();
    for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
      unsigned type = graph.bond_type(b);
      unsigned site_i = graph.source(b);
      unsigned site_j = graph.target(b);
      int bc_phase = graph.bond_sign(b);
      // matrix elements each term & bond type
      unsigned term = 0;
      for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
        matrix_elem(term,type) += config.apply(it->qn_operator(),site_i,site_j,bc_phase);
        term++;
      }
    }
    unsigned i = 0;
    for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
      for (unsigned btype=0; btype<graph.num_bond_types(); ++btype) {
        config_energy_(i) += std::real(it->coupling(btype)*matrix_elem(i,btype));
      }
      i++;
    }
  }

  // site energies
  if (model.has_siteterm()) {
    Matrix matrix_elem(model.num_siteterms(),graph.num_site_types());
    matrix_elem.setZero();
    Eigen::VectorXi hubbard_nd(graph.num_site_types());
    hubbard_nd.setZero();
    // do it little differently
    unsigned term = model.num_bondterms();
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
          unsigned site = graph.site(s);
          unsigned type = graph.site_type(s);
          hubbard_nd(type) += config.apply_niup_nidn(site);
        }
      }
      else {
        for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
          unsigned site = graph.site(s);
          unsigned type = graph.site_type(s);
          matrix_elem(term,type) += config.apply(it->qn_operator(), site);
        }
      }
      term++;
    }
    /*for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
      unsigned type = graph.site_type(s);
      // matrix elements each term & site type
      unsigned term = model.num_bondterms();
      for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
        //matrix_elem(term,type) += config.apply(it->qn_operator(), s);
        term++;
      }
    }*/
    unsigned i = 0;
    unsigned n = model.num_bondterms();
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          config_energy_(n+i) += std::real(it->coupling(stype)*hubbard_nd(stype));
        }
      }
      else {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          config_energy_(n+i) += std::real(it->coupling(stype)*matrix_elem(i,stype));
        }
      }
      i++;
    }
  }
  // energy per site
  return config_energy_/num_sites_;
}



} // end namespace vmc

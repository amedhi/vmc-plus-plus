/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-17 23:30:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-03 23:00:31
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include "simulator.h"

namespace vmc {

int Simulator::do_measurements(void)
{
  if (optimizing_mode_) {
    double energy = config_energy().sum();
    observables.total_energy() << energy;
    //std::cout << "--------here--------\n";
    if (observables.energy_grad()) {
      config.get_grad_logpsi(grad_logpsi_);
      unsigned n = 0;
      for (unsigned i=0; i<num_varparms_; ++i) {
        energy_grad2_[n] = energy * grad_logpsi_[i];
        energy_grad2_[n+1] = grad_logpsi_[i];
        n += 2;
      }
      observables.energy_grad2() << energy_grad2_;
    }
    if (observables.sr_matrix()) {
      if (!observables.energy_grad()) 
        throw std::logic_error("Simulator::do_measurements: internal error");
      // operator 'del(ln(psi))' terms
      for (unsigned i=0; i<num_varparms_; ++i) 
        sr_matrix_el_[i] = grad_logpsi_[i];
      // flatten the upper triangular part to a vector
      unsigned k = num_varparms_;
      for (unsigned i=0; i<num_varparms_; ++i) {
        double x = grad_logpsi_[i];
        for (unsigned j=i; j<num_varparms_; ++j) {
          double y = grad_logpsi_[j];
          sr_matrix_el_[k] = x * y;
          ++k;
        }
      }
      observables.sr_matrix() << sr_matrix_el_;
    }
    return 0;
  }

  // normal run
  if (need_energy_) {
    term_energy_ = config_energy();
    //std::cout << "--------here--------\n";
    if (observables.energy()) observables.energy() << term_energy_;

    if (observables.energy_grad()) {
      double energy = term_energy_.sum();
      observables.total_energy() << energy;
      config.get_grad_logpsi(grad_logpsi_);
  //-------------------------------------------------
      unsigned n = 0;
      for (unsigned i=0; i<num_varparms_; ++i) {
        energy_grad2_[n] = energy * grad_logpsi_[i];
        energy_grad2_[n+1] = grad_logpsi_[i];
        n += 2;
      }
      observables.energy_grad2() << energy_grad2_;
    } 

    else if (observables.total_energy()) {
      double energy = term_energy_.sum();
      observables.total_energy() << energy;
    }
  }
  /*if (observables_.total_energy()) {
    double e = config_energy().sum();
    observables_.total_energy() << e;
    observables.total_energy() << e;
  }*/
  return 0;
}

int Simulator::finalize_energy_grad(void)
{
  double mean_energy = observables.total_energy().mean();
  energy_grad2_ = observables.energy_grad2().mean_data(); 
  unsigned n = 0;
  for (unsigned i=0; i<num_varparms_; ++i) {
    energy_grad_(i) = 2.0*(energy_grad2_[n] - mean_energy * energy_grad2_[n+1]);
    n += 2;
  }
  //-------------------------------------------------
  //std::cout << mean_energy << "\n";
  //std::cout << energy_grad2_ << "\n";
  //for (int i=0; i<energy_grad_.size(); ++i) std::cout << energy_grad_[i] << "\n";
  //-------------------------------------------------
  observables.energy_grad() << energy_grad_;
  return 0;
}


obs::vector Simulator::config_energy(void) const
{
  using op_id = model::op_id;
  //for (auto& elem : term_energy_) elem = 0.0;
  term_energy_.setZero();
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
        term_energy_(i) += std::real(it->coupling(btype)*matrix_elem(i,btype));
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
          term_energy_(n+i) += std::real(it->coupling(stype)*hubbard_nd(stype));
        }
      }
      else {
        for (unsigned stype=0; stype<graph.num_site_types(); ++stype) {
          term_energy_(n+i) += std::real(it->coupling(stype)*matrix_elem(i,stype));
        }
      }
      i++;
    }
  }
  // energy per site
  return term_energy_/num_sites_;
}



} // end namespace vmc

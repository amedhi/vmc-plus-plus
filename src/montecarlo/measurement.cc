/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include <set>
#include <cmath>
#include "simulator.h"

namespace mc {

/*----------------------Energy-----------------------*/
mc::VectorData Simulator::get_energy(void)
{
  unsigned num_terms = num_total_terms();
  mc::VectorData state_energy(num_terms);
  for (unsigned i=0; i<num_terms; ++i) state_energy(i) = 0.0;
  // bond energies
  for (auto b=bonds_begin(); b!=bonds_end(); ++b) {
    unsigned type = bond_type(b);
    auto src = source(b);
    auto tgt = target(b);
    state_idx src_idx = state[site(src)].idx();
    state_idx tgt_idx = state[site(tgt)].idx();
    unsigned term = 0;
    for (auto bterm=bondterms_begin(); bterm!=bondterms_end(); ++bterm) {
      double m = bterm->matrix_element(type, src_idx, tgt_idx);
      double c = bterm->coupling(type);
      //std::cout << "bond: " << " " << m << std::endl;
      state_energy[term++] += m * c;
    }
  }
  // site energies
  for (const auto& s : state) {
    unsigned term = num_bondterms();
    for (auto sterm=siteterms_begin(); sterm!=siteterms_end(); ++sterm) {
      double m = sterm->matrix_element(s.type(), s.idx());
      double c = sterm->coupling(s.type());
      state_energy[term++] += m * c;
    }
  }
  // energy per site
  return state_energy/num_sites();
}

/*----------------------Ising magnetization-----------------------*/
double Simulator::get_magnetization(void)
{
  double ms = 0;
  for (const auto& s : state) {
    ms += observables.magn_op().apply(s);
  }
  return std::abs(ms)/num_sites();
};

/*----------------------Potts magnetization-----------------------*/
double Simulator::get_potts_magnetization(void)
{
  unsigned num_sitetypes = basis().size();
  std::vector<std::multiset<int> > lattice_spins(num_sitetypes);
  std::vector<std::set<int> > site_spins(num_sitetypes);

  // spin values for a sites
  for (unsigned i=0; i<num_sitetypes; ++i) {
    for (unsigned j=0; j<sitebasis_dimension(i); ++j) {
      SiteBasisState state(i, j, sitebasis_dimension(i)-1); 
      int spin = std::nearbyint(observables.potts_magn_op().apply(state));
      site_spins[i].insert(spin);
    }
  }

  // set of spin values for all the sites in the lattice
  for (const auto& s : state) {
    int spin = std::nearbyint(observables.potts_magn_op().apply(s));
    lattice_spins[s.type()].insert(spin);
  }

  // magnetization
  double ms = 0;
  for (unsigned i=0; i<num_sitetypes; ++i) {
    int q = site_spins[i].size();
    int Nmax = 0;
    for (const auto& s : site_spins[i]) {
      int count = lattice_spins[i].count(s);
      if (Nmax < count) Nmax = count;
    }
    ms += static_cast<double>(q*Nmax - lattice_spins[i].size())/(q-1);
  }
  return ms/num_sites();
}

/*----------------------Strain-----------------------*/
double Simulator::get_strain(void)
{
  double ms = 0;
  for (const auto& s : state) {
    ms += observables.strain_op().apply(s);
  }
  return std::abs(ms)/num_sites();
};

} // end namespace basis

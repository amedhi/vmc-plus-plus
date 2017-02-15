/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:20:28
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-15 22:41:54
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "basisstate.h"
#include <stdexcept>
#include <algorithm>

namespace vmc {

/*
unsigned Hole::total_num_ = 0;
unsigned SpinUp::total_num_ = 0;
unsigned SpinDn::total_num_ = 0;
unsigned Doublon::total_num_ = 0;
*/

BasisState::BasisState(void) 
{
  set_vaccuum(0);
} 

BasisState::BasisState(const unsigned& num_sites) 
{
  set_vaccuum(num_sites);
} 

BasisState::BasisState(const unsigned& num_sites, const bool& allow_dbl) 
{
  set_vaccuum(num_sites, allow_dbl);
} 

void BasisState::set_vaccuum(const unsigned& num_sites, const bool& allow_dbl) 
{
  num_sites_ = num_sites;
  resize(num_sites_);
  double_occupancy_ = allow_dbl;
  for (unsigned i=0; i<num_sites_; ++i) {
    operator[](i).put_uphole(i);
    operator[](i).put_dnhole(i);
  }
  upspin_pos_.clear();
  dnspin_pos_.clear();
}

void BasisState::clear(void)
{
  upspin_pos_.clear();
  dnspin_pos_.clear();
} 

void BasisState::allow_double_occupancy(const bool& allow)
{
  double_occupancy_ = allow;
} 

bool BasisState::create_upspin(const unsigned& site)
{
  //this->at(site).put_upspin(upspins_.size());
  //upspins_.push_back(&at(site));
  return true;
}

void BasisState::init_random(const unsigned& num_upspins, const unsigned& num_dnspins,
    RandomGenerator& rng, const bool& allow_dbl)
{
  if (num_upspins>num_sites_ || num_dnspins>num_sites_)
    throw std::range_error("* BasisState::create_random: spin number exceed capacity");
  upspin_pos_.resize(num_upspins);
  dnspin_pos_.resize(num_dnspins);
  uphole_pos_.resize(num_sites_ - num_upspins);
  dnhole_pos_.resize(num_sites_ - num_dnspins);
  double_occupancy_ = allow_dbl;

  std::vector<unsigned> all_sites(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) all_sites[i] = i;
  std::shuffle(all_sites.begin(),all_sites.end(),rng);

  // UP spins & holes
  for (unsigned i=0; i<num_upspins; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_upspin(i);
    upspin_pos_[i] = site;
  }
  unsigned uh = 0;
  for (unsigned i=num_upspins; i<num_sites_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_uphole(uh);
    uphole_pos_[uh] = site;
    uh++;
  }
  // DN spins & holes
  if (double_occupancy_) {
    std::shuffle(all_sites.begin(),all_sites.end(),rng);
    for (unsigned i=0; i<num_dnspins; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnspin(i);
      dnspin_pos_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=num_dnspins; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(dh);
      dnhole_pos_[dh] = site;
      dh++;
    }
  }
  else {
    unsigned total_spins = num_upspins+num_dnspins;
    if (total_spins > num_sites_)
      throw std::range_error("* BasisState::create_random: spin number exceed capacity");
    // DN spins
    unsigned ds = 0;
    for (unsigned i=num_upspins; i<total_spins; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnspin(ds);
      dnspin_pos_[ds] = site;
      ds++;
    }
    // DN holes
    for (unsigned i=0; i<num_upspins; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_uphole(i);
      uphole_pos_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=total_spins; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(num_upspins+dh);
      dnhole_pos_[num_upspins+dh] = site;
      dh++;
    }
  }
}

void BasisState::hop_upspin(const unsigned& i, const unsigned& s)
{
  //at(upspin_pos[i]).reset(spin::UP);
  //at(s).create_upspin(i);
  //upspin_pos[i] = s;
}

std::ostream& operator<<(std::ostream& os, const BasisState& bs)
{
  unsigned len = 12;
  unsigned i = 0;
  os << std::string(60,'-') << std::endl;
  os << "Basis state:\n";
  for (const auto& s : bs) {
    os << "(" << s.to_string() << ") ";
    i++;
    if (i==len) {
      os << std::endl; i = 0;
    }
  }
  os << std::endl; 
  os << std::string(60,'-') << std::endl;
  return os;
}



} // end namespace vmc

/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:20:28
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 13:04:26
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
  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
}

void BasisState::clear(void)
{
  upspin_pos_.clear();
  dnspin_pos_.clear();
} 

void BasisState::init_spins(const unsigned& num_upspins, const unsigned& num_dnspins, 
    const bool& allow_dbloccupancy)
{
  num_upspins_ = num_upspins;
  num_dnspins_ = num_dnspins;
  double_occupancy_ = allow_dbloccupancy;
  if (num_upspins_>num_sites_ || num_dnspins_>num_sites_)
    throw std::range_error("* BasisState::init_spins: spin number exceed capacity");
  if (!double_occupancy_ && (num_upspins_+num_dnspins_)>num_sites_)
    throw std::range_error("* BasisState::init_spins: spin number exceed capacity");
  num_upholes_ = num_sites_ - num_upspins;
  num_dnholes_ = num_sites_ - num_dnspins;
  // resizing
  upspin_pos_.resize(num_upspins_);
  dnspin_pos_.resize(num_dnspins_);
  uphole_pos_.resize(num_upholes_);
  dnhole_pos_.resize(num_dnholes_);
  // random generator
  rng_.set_upspin_generator(0,num_upspins_-1);
  rng_.set_dnspin_generator(0,num_dnspins_-1);
  rng_.set_uphole_generator(0,num_upholes_-1);
  rng_.set_dnhole_generator(0,num_dnholes_-1);
} 

void BasisState::allow_double_occupancy(const bool& allow)
{
  double_occupancy_ = allow;
} 

void BasisState::set_random(void)
{
  std::vector<unsigned> all_sites(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) all_sites[i] = i;
  std::shuffle(all_sites.begin(),all_sites.end(),rng_);

  // UP spins & holes
  for (unsigned i=0; i<num_upspins_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_upspin(i);
    upspin_pos_[i] = site;
  }
  unsigned uh = 0;
  for (unsigned i=num_upspins_; i<num_sites_; ++i) {
    unsigned site = all_sites[i];
    operator[](site).put_uphole(uh);
    uphole_pos_[uh] = site;
    uh++;
  }
  // DN spins & holes
  if (double_occupancy_) {
    std::shuffle(all_sites.begin(),all_sites.end(),rng_);
    for (unsigned i=0; i<num_dnspins_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnspin(i);
      dnspin_pos_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=num_dnspins_; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(dh);
      dnhole_pos_[dh] = site;
      dh++;
    }
  }
  else {
    unsigned total_spins = num_upspins_+num_dnspins_;
    // DN spins
    unsigned ds = 0;
    for (unsigned i=num_upspins_; i<total_spins; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnspin(ds);
      dnspin_pos_[ds] = site;
      ds++;
    }
    // DN holes
    for (unsigned i=0; i<num_upspins_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_uphole(i);
      uphole_pos_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=total_spins; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(num_upspins_+dh);
      dnhole_pos_[num_upspins_+dh] = site;
      dh++;
    }
  }
}

std::pair<int,int> BasisState::random_upspin_hop(void)
{
  //at(upspin_pos[i]).reset(spin::UP);
  //at(s).create_upspin(i);
  //upspin_pos[i] = s;
  return std::make_pair(0,0);
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

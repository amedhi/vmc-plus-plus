/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:20:28
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-16 17:17:57
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
  uphole_pos_.clear();
  dnhole_pos_.clear();
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
  // assuming 'num_upspins>0', 'num_dnspins>0'
  rng_.set_upspin_generator(0,static_cast<int>(num_upspins_)-1);
  rng_.set_dnspin_generator(0,static_cast<int>(num_dnspins_)-1);
  // hole numbers may be zero
  int m = std::max(static_cast<int>(num_upholes_),1);
  rng_.set_uphole_generator(0,m-1);
  int n = std::max(static_cast<int>(num_dnholes_),1);
  if (num_dnholes_>0)
    rng_.set_dnhole_generator(0,n-1);
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
      operator[](site).put_dnhole(i);
      dnhole_pos_[i] = site;
    }
    unsigned dh = 0;
    for (unsigned i=total_spins; i<num_sites_; ++i) {
      unsigned site = all_sites[i];
      operator[](site).put_dnhole(num_upspins_+dh);
      dnhole_pos_[num_upspins_+dh] = site;
      dh++;
    }
  }
  // in addition, in case of 'no double occupancy':
  if (!double_occupancy_) {
    for (unsigned i=0; i<num_upspins_; ++i) {
      if (operator[](uphole_pos_[i]).have_dnspin())
        uphole_pos_[i] = -1;
    }
    for (unsigned i=0; i<num_dnspins_; ++i) {
      if (operator[](dnhole_pos_[i]).have_upspin())
        dnhole_pos_[i] = -1;
    }
  }
  // in case there are no holes
  if (num_upholes_==0) {
    uphole_pos_.resize(1);
    upspin_pos_[0] = -1;
  }
  if (num_dnholes_==0) {
    dnhole_pos_.resize(1);
    dnspin_pos_[0] = -1;
  }
}

const std::pair<int,int>& BasisState::random_upspin_hop(void)
{
  spin_site_pair.first = rng_.random_upspin();
  move_hole = rng_.random_uphole();
  spin_site_pair.second = uphole_pos_[move_hole]; 
  if (spin_site_pair.second < 0) 
    proposed_move=move_type::null;
  else 
    proposed_move=move_type::upspin_hop;
  /*if (double_occupancy_) {
    proposed_move = move_type::upspin_hop;
  }
  else {
    proposed_move = move_type::null;
    spin_site_pair.second = -1;
  }
  */
  return spin_site_pair;
}

void BasisState::accept_last_move(void)
{
  int site;
  switch (proposed_move) {
    case move_type::upspin_hop:
      // hole moves back
      site = upspin_pos_[spin_site_pair.first];
      operator[](site).put_uphole(move_hole);
      uphole_pos_[move_hole] = site;
      // remove block for dnhole (in case no double occupancy)
      if (operator[](site).have_dnhole())
        dnhole_pos_[operator[](site).dnhole_id()] = site;
      // spin moves forward
      operator[](spin_site_pair.second).put_upspin(spin_site_pair.first);
      upspin_pos_[spin_site_pair.first] = spin_site_pair.second;
      proposed_move = move_type::null;
      break;
    case move_type::dnspin_hop:
      // hole moves back
      site = dnspin_pos_[spin_site_pair.first];
      operator[](site).put_dnhole(move_hole);
      dnhole_pos_[move_hole] = site;
      // remove block for dnhole (in case no double occupancy)
      if (operator[](site).have_uphole())
        uphole_pos_[operator[](site).uphole_id()] = site;
      // spin moves forward
      operator[](spin_site_pair.second).put_dnspin(spin_site_pair.first);
      dnspin_pos_[spin_site_pair.first] = spin_site_pair.second;
      proposed_move = move_type::null;
      break;
    case move_type::exchange:
      break;
    case move_type::null:
      break;
  }
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

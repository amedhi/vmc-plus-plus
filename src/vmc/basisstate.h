/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:16:02
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-15 23:01:38
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef BASISSTATE_H
#define BASISSTATE_H

#include <iostream>
#include <vector>
#include <utility>
#include <bitset>
#include <array>
#include "./random.h"

namespace vmc {

enum {MAX_SITES=1000};

enum class spin {UP, DN};

class SiteState : public std::bitset<2>
{
public:
  SiteState()
    : null_id_{100000}  
    , id_{null_id_,null_id_}
    , spin_UP_{static_cast<size_t>(spin::UP)}
    , spin_DN_{static_cast<size_t>(spin::DN)}
  { 
    reset(); 
  }
  void put_upspin(const unsigned& n) {set(spin_UP_); id_.first=n;} 
  void put_uphole(const unsigned& n) {reset(spin_UP_); id_.first=n;} 
  void put_dnspin(const unsigned& n) {set(spin_DN_); id_.second=n;} 
  void put_dnhole(const unsigned& n) {reset(spin_DN_); id_.second=n;} 
  const unsigned& upspin_id(void) const { return id_.first; }
  const unsigned& dnspin_id(void) const { return id_.second; }
  bool have_upspin(void) const { return test(spin_UP_); }  
  bool have_dnspin(void) const { return test(spin_DN_); }  
  bool is_full(void) const { return all(); }  
  bool is_empty(void) const { return none(); }  
  bool is_not_empty(void) const { return any(); }  
private:
  using size_t = std::size_t;
  unsigned null_id_{100000};
  std::pair<unsigned,unsigned> id_;
  unsigned spin_UP_;
  unsigned spin_DN_;
};

class BasisState : public std::vector<SiteState> 
{
public:
  BasisState(); 
  BasisState(const unsigned& num_sites); 
  BasisState(const unsigned& num_sites, const bool& allow_dbl);
  ~BasisState() {} 
  void set_vaccuum(const unsigned& num_sites, const bool& allow_dbl=true);
  void allow_double_occupancy(const bool& allow);
  void init_random(const unsigned& num_upspins, const unsigned& num_dnspins,
    RandomGenerator& rng, const bool& allow_double_occupancy=true);
  bool create_upspin(const unsigned& site);
  bool create_dnspin(const unsigned& site);
  void init(const unsigned& num_sites); 
  void hop_upspin(const unsigned& i, const unsigned& s);
  const std::vector<unsigned>& upspin_pos(void) const { return upspin_pos_; }
  const std::vector<unsigned>& dnspin_pos(void) const { return dnspin_pos_; }
  //const SiteState& upspin(unsigned)
  //bool annihilate(i,sigma);
  //bool hop(i,j,sigma);
  friend std::ostream& operator<<(std::ostream& os, const BasisState& bs);
private:
  enum {max_sites=1000};
  unsigned num_sites_{0};
  bool double_occupancy_{true};
  std::vector<unsigned> upspin_pos_;
  std::vector<unsigned> dnspin_pos_;
  std::vector<unsigned> uphole_pos_;
  std::vector<unsigned> dnhole_pos_;
  void clear(void); 
};


/*
class SpinParticle 
{
public:
  SpinParticle() : id_{0}, site_{0}, sigma_{spin::HL} {}
  SpinParticle(const spin& s) : id_{0}, site_{0}, sigma_{s} {}
  SpinParticle(const unsigned& i, const spin& s) : id_{0}, site_{i}, sigma_{s} {}
  const spin& sigma(void) const { return sigma_; }
  const unsigned& id(void) const { return id_; }
protected:
  unsigned id_;
private:
  unsigned site_;
  spin sigma_;
};

class SpinUp : public SpinParticle 
{
public:
  SpinUp() {}
  SpinUp(const unsigned& site) : SpinParticle{site, spin::UP} { id_=total_num_++; }
  SpinUp(const SpinUp& p) : SpinParticle(p) { total_num_++; }
  ~SpinUp() { --total_num_; }
private:
  static unsigned total_num_;
};

class SpinDn : public SpinParticle 
{
public:
  SpinDn(const unsigned& site) : SpinParticle{site, spin::DN} { id_=total_num_++; }
  SpinDn(const SpinDn& p) : SpinParticle(p) { total_num_++; }
  ~SpinDn() { --total_num_; }
  static unsigned total_num_;
private:
};

class Hole : public SpinParticle 
{
public:
  Hole(const unsigned& site) : SpinParticle{site, spin::HL} { id_=total_num_++; }
  Hole(const Hole& p) : SpinParticle(p) { total_num_++; }
  ~Hole() { --total_num_; }
  static unsigned num_holes(void) { return total_num_; }
private:
  static unsigned total_num_;
};

class Doublon : public SpinParticle 
{
public:
  Doublon() : SpinParticle(spin::UD) { id_=total_num_++; }
  ~Doublon() { --total_num_; }
private:
  unsigned up_id_;
  unsigned dn_id_;
  static unsigned total_num_;
};
using site_state = SpinParticle*;

*/




} // end namespace vmc

#endif
/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef RANDOM_H
#define RANDOM_H

#include <iostream>
#include <stdexcept>
#include <chrono>
#include <array>
#include <map>
#include <random>
#include <time.h> 
#include "../lattice/lattice.h" 

namespace mc {

class RandomNumber : public std::mt19937_64
{
public:
  RandomNumber();
  RandomNumber(const unsigned& seed_type);
  ~RandomNumber() {};
  void set_site_generator(const unsigned& min, const unsigned& max);
  void add_state_generator(const unsigned& site_type, const unsigned& min, const unsigned& max);
  void seed(const int& seed_type);
  void time_seed(void);
  //unsigned random_idx(const unsigned& site_type) {return state_dist_map[site_type](*this); }
  unsigned random_idx(const unsigned& site_type) { return state_generators[site_type](*this); }
  //unsigned random_site(void) { return site_dist[0](*this); }
  unsigned random_site(void) { return site_generator(*this); }
  double random_real(void) { return real_generator(*this); }
private:
  using int_dist = std::uniform_int_distribution<unsigned>;
  using myclock = std::chrono::high_resolution_clock;
  int seed_type_;
  //std::map<unsigned, int_dist> site_dist; // there will be only one item here
  int_dist site_generator; // there will be only one item here
  //std::map<unsigned, int_dist> state_dist_map;
  std::array<int_dist, lattice::MAX_SITE_TYPES> state_generators;
  std::uniform_real_distribution<double> real_generator;
};


} // end namespace monte carlo

#endif

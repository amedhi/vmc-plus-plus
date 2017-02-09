/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2016-03-16 16:47:10
*----------------------------------------------------------------------------*/
#include "random.h"

namespace mc {

RandomNumber::RandomNumber() : seed_type_(0), real_generator(0.0, 1.0) 
{
  for (auto& g : state_generators) g = int_dist(-1,-1);
}
  
RandomNumber::RandomNumber(const unsigned& seed_type) 
  : seed_type_(seed_type), real_generator(0.0, 1.0)
{
  if (seed_type_==1) time_seed();
  for (auto& g : state_generators) g = int_dist(-1,-1);
}

void RandomNumber::seed(const int& seed_type) {
  seed_type_ = seed_type;
  if (seed_type_==1) time_seed();
} 

void RandomNumber::time_seed(void) 
{
  myclock::time_point now = myclock::now();
  myclock::duration till_now = now.time_since_epoch();
  unsigned itc = till_now.count();
  this->std::mt19937_64::seed(itc);
}

void RandomNumber::add_state_generator(const unsigned& site_type, const unsigned& min, 
  const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomNumber::add_state_generator: invalid input");
  state_generators[site_type] = int_dist(min, max);
  /*if (state_dist_map.find(site_type) != state_dist_map.end()) return;
  // nothing to be done, the dist already exists. Else:
  if (min>max) throw std::runtime_error("RandomNumber::add_state_generator: invalid input");
  state_dist_map.insert(std::make_pair(site_type, int_dist(min, max)));
  return;
  */
}

void RandomNumber::set_site_generator(const unsigned& min, const unsigned& max)
{
  if (min>max) throw std::runtime_error("RandomNumber::set_site_generator: invalid input");
  site_generator = int_dist(min, max);
  /*if (site_dist.size() == 0) 
    site_dist.insert(std::make_pair(0, int_dist(min, max)));
  else
    site_dist[0] = int_dist(min, max);
  */
}








} // end namespace basis

/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef SITEBASISSTATE_H
#define SITEBASISSTATE_H


namespace mc {

class SiteBasisState
{
public:
  using state_idx = unsigned;
  SiteBasisState() : type_{0}, idx_{0}, max_idx_{0} {}
  SiteBasisState(const unsigned& type, const state_idx& idx, const state_idx& max_idx) 
    : type_{type}, idx_{idx}, max_idx_{max_idx} {}
  SiteBasisState(const unsigned& type, const state_idx& max_idx) 
    : type_{type}, idx_{0}, max_idx_{max_idx} {}
  ~SiteBasisState() {};

  inline const state_idx& type(void) const { return type_; }
  inline const state_idx& idx(void) const { return idx_; }
  inline const state_idx& max_idx(void) const { return max_idx_; }
  inline void reset_idx(const state_idx& idx) { idx_=idx; }
private:
  unsigned type_{0};
  state_idx idx_{0};
  state_idx max_idx_{0}; 
};

} // end namespace monte carlo

#endif

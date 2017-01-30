/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef OBSERVABLE_OPERATOR_H
#define OBSERVABLE_OPERATOR_H

#include <string>
#include <vector>
#include <set>
#include <stdexcept>
#include <Eigen/Core>
#include "../model/hamiltonian_term.h"
#include "../model/model.h"
#include "../montecarlo/sitebasisstate.h"

namespace mc {

class SiteObsOperator : public std::vector<model::SiteOperator> 
{
  using op_sitetypes = std::set<unsigned>;
public:
  SiteObsOperator() {}
  SiteObsOperator(const model::BasisDescriptor& basis, const std::string& op_expr, 
    const std::string& site="i");
  SiteObsOperator(const model::BasisDescriptor& basis, const op_sitetypes& sitetypes, 
    const std::string& op_expr, const std::string& site="i");
  ~SiteObsOperator() {}
  void init(const model::BasisDescriptor& basis, const std::string& op_expr, 
    const std::string& site="i");
  void init(const model::BasisDescriptor& basis, const op_sitetypes& sitetypes, 
    const std::string& op_expr, const std::string& site="i");
  inline const double& apply(const SiteBasisState& state) const
  { return operator[](state.type()).matrix_element(state.idx()); }
private:
};


} // end namespace mc

#endif

/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observable_operator.h"
#include "../expression/expression.h"

namespace mc {

SiteObsOperator::SiteObsOperator(const model::BasisDescriptor& basis, const std::string& op_expr, 
    const std::string& site)
{
  init(basis, op_expr, site);
}

SiteObsOperator::SiteObsOperator(const model::BasisDescriptor& basis, 
  const op_sitetypes& sitetypes, const std::string& op_expr, const std::string& site)
{
  init(basis, sitetypes, op_expr, site);
}

void SiteObsOperator::init(const model::BasisDescriptor& basis, const std::string& op_expr, 
    const std::string& site)
{
  // sitetype not specied, acts on all sites
  // set all the symbols in the 'op_expr' to default value zero
  expr::Expression::variables vars;
  // strip '(' & ')' characters around 'site' in op_expr
  std::string expr_str(op_expr);
  std::string s = "("+site+")";
  std::string::size_type pos;
  while ((pos=expr_str.find(s)) != std::string::npos) expr_str.replace(pos,s.length(),site);

  for (const auto& elem : basis) {
    unsigned num_operator = elem.second.num_operator();
    for (unsigned n=0; n<num_operator; ++n) {
      model::QuantumOperator op = elem.second.quantum_operator(n);
      std::string vname = op.name()+site;
      vars[vname] = 0.0;
    }
  }

  // now build the matrices
  for (unsigned type=0; type<basis.size(); ++type) {
    push_back(model::SiteOperator(op_expr,site));
    back().build_matrix(basis.at(type), vars);
  }

}

void SiteObsOperator::init(const model::BasisDescriptor& basis, 
  const op_sitetypes& sitetypes, const std::string& op_expr, const std::string& site)
{
  // construct the matrices for specified site_types
  for (unsigned type=0; type<basis.size(); ++type) {
    if (sitetypes.find(type) != sitetypes.end()) {
      push_back(model::SiteOperator(op_expr,site));
    }
    else push_back(model::SiteOperator());
    back().build_matrix(basis.at(type));
  }
}

} // end namespace mc


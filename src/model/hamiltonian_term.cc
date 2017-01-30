/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   amedhi
* Last Modified time: 2017-01-30 19:12:33
*----------------------------------------------------------------------------*/
#include "hamiltonian_term.h"

namespace model {

/*--------------------------CouplingConstant--------------------*/
const int CouplingConstant::global_type = -1;

CouplingConstant::CouplingConstant(const std::string& expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr});
  num_types_ = -1;
  valid_ = true;
}

std::pair<CouplingConstant::iterator, bool> CouplingConstant::insert(const value_type& val) 
{
  std::pair<iterator, bool> res = super_type::insert(val);
  num_types_ = super_type::size();
  valid_ = true;
  return res;
}

CouplingConstant& CouplingConstant::operator=(const std::string expr)
{
  super_type::clear();
  // expr is applicable for all site & bond types
  super_type::insert({global_type, expr});
  num_types_ = -1; 
  valid_ = true;
  return *this;
}

CouplingConstant::CouplingConstant(const value_type& type0, const value_type& type1, 
    const value_type& type2, const value_type& type3, const value_type& type4, 
    const value_type& type5)
{
  create(type0, type1, type2, type3, type4, type5); 
}

void CouplingConstant::clear(void)
{
  super_type::clear();
  num_types_ = 0;
  valid_ = false;
}

void CouplingConstant::create(const unsigned& num_types) 
{
  super_type::clear();
  num_types_ = num_types;
  valid_ = false;
}

void CouplingConstant::create(const value_type& type0, const value_type& type1, 
    const value_type& type2, const value_type& type3, const value_type& type4, 
    const value_type& type5)
{
  super_type::insert(type0);
  num_types_=1;
  if (type1.second != "_null_") {
    super_type::insert(type1); num_types_++;
  }
  if (type2.second != "_null_") {
    super_type::insert(type2); num_types_++;
  }
  if (type3.second != "_null_") {
    super_type::insert(type3); num_types_++;
  }
  if (type4.second != "_null_") {
    super_type::insert(type4); num_types_++;
  }
  if (type5.second != "_null_") {
    super_type::insert(type5); num_types_++;
  }
  valid_=true;
}

void CouplingConstant::add_type(const unsigned& type, const std::string& expr)
{
  super_type::insert({type, expr});
  valid_ = (num_types_==static_cast<int>(size()));
}

void CouplingConstant::add_type(const value_type& val)
{
  super_type::insert(val);
  valid_ = (num_types_==static_cast<int>(size()));
}


//-----------------------SiteOperatorTerm-------------------------
SiteOperatorTerm::SiteOperatorTerm(const std::string& name, const std::string& cc_expr, 
  const qn_op& op)
  : op_{op}, name_{name}, cc_expr_{cc_expr}
{
}

int SiteOperatorTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  if (cc_expr_.size()==0) {
    cc_value_ = 0.0;
    return 0;
  }

  // expression evaluator
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& c : cvals) vars[c.first] = c.second;
  for (const auto& p : pvals) vars[p.first] = p.second;
  try { 
    cc_value_ = expr.evaluate(cc_expr_, vars); 
  }
  catch (std::exception& e) 
  { 
    std::string msg = "SiteOperatorTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
  return 0;
}

SiteTerm::SiteTerm(const std::string& name, const CouplingConstant& cc, 
    const qn_op& op, const unsigned& size)
{
  if (!cc.valid()) throw std::invalid_argument("SiteTerm::Invalid CouplingConstant");
  name_ = name;
  size_ = size;

  // if the 'cc' is implicitly defined for all site types 
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    std::string cc_expr = cc.begin()->second;
    for (unsigned i=0; i<size_; ++i) {
      std::string term_name = name + std::to_string(i);
      this->operator[](i) = SiteOperatorTerm(term_name,cc_expr,op);
    }
  } 
  else {
    // set "SiteOperatorTerm"-s only for those site types for which 'cc' is 
    // explicitly defined
    for (const auto& p : cc) {
      std::string term_name = name + std::to_string(p.first);
      std::string cc_expr = p.second;
      //std::cout << p.first << " " << p.second << "\n";
      this->operator[](p.first) = SiteOperatorTerm(term_name,cc_expr,op);
      //insert({p.first, SiteOperatorTerm(term_name,p.second,op_expr,site)});
    }
  }
}

void SiteTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  for (auto& term : *this) term.eval_coupling_constant(cvals, pvals);
}

const double& SiteTerm::coupling(const unsigned& site_type) const
{
  return this->operator[](site_type).coupling();
}


//-----------------------BondOperatorTerm-------------------------
BondOperatorTerm::BondOperatorTerm(const std::string& name, const std::string& cc_expr, 
  const qn_op& op)
  : op_{op}, name_{name}, cc_expr_{cc_expr}, cc_value_{0.0}
{
}

int BondOperatorTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  if (cc_expr_.size()==0) {
    cc_value_ = 0.0;
    return 0;
  }
  // expression evaluator
  expr::Expression expr;
  expr::Expression::variables vars;
  for (const auto& c : cvals) vars[c.first] = c.second;
  for (const auto& p : pvals) vars[p.first] = p.second;
  try { 
    cc_value_ = expr.evaluate(cc_expr_, vars); 
    //std::cout << "bterm: " << name_ << " = " << cc_value_ << "\n";
  }
  catch (std::exception& e) 
  { 
    std::string msg = "BondOperatorTerm::evaluate_coupling_constant:\n" + std::string(e.what());
    throw std::runtime_error(msg);
  }
  return 0;
}

//-----------------------BondTerm-------------------------
BondTerm::BondTerm(const std::string& name, const CouplingConstant& cc, const qn_op& op,
    const unsigned& size)
{
  if (!cc.valid()) throw std::invalid_argument("BondTerm:: Invalid CouplingConstant");
  name_ = name;
  size_ = size;

  // if the 'cc' is implicitly defined for all bond types 
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    std::string cc_expr = cc.begin()->second;
    for (unsigned i=0; i<size_; ++i) {
      std::string term_name = name + std::to_string(i);
      this->operator[](i) = BondOperatorTerm(term_name,cc_expr,op);
    }
  } 
  else {
    // set "BondOperatorTerm"-s only for those bond types for which 'cc' is 
    // defined explicitly
    for (const auto& p : cc) {
      std::string term_name = name + std::to_string(p.first);
      std::string cc_expr = p.second;
      this->operator[](p.first) = BondOperatorTerm(term_name,cc_expr,op);
      //insert({p.first, BondOperatorTerm(term_name,p.second,op_expr,src,tgt)});
    }
  }
}

void BondTerm::eval_coupling_constant(const ModelParams& cvals, const ModelParams& pvals)
{
  for (auto& term : *this) term.eval_coupling_constant(cvals, pvals);
  //std::cout << "------hi--------\n"; //abort();
  /*for (auto& elem : *this) {
    elem.second.eval_coupling_constant(pvals);
  }*/
}

const double& BondTerm::coupling(const unsigned& bond_type) const
{
  return this->operator[](bond_type).coupling();
  /*const_iterator it=find(bond_type);
  if (it != end()) return it->second.coupling();
  else return null_coupling_;
  */
}



} // end namespace model

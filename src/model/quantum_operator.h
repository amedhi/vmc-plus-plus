/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef QUANTUM_OPERATOR_H
#define QUANTUM_OPERATOR_H

#include <string>
#include <vector>
#include "qn.h"

namespace model {


class QuantumOperator : private std::vector<std::pair<std::string,QuantumNumber::value_type> >
{
public:
  using QN = QuantumNumber;
  //QuantumOperator();
  QuantumOperator(const std::string& name, const std::string& operand_name);
  QuantumOperator(const std::string& name, const std::string& elem, 
    const std::string& operand_name, const QN::value_type& change);
  ~QuantumOperator() {};
  unsigned add_operand(const std::string& operand_name, const short& change) 
  { 
    this->push_back(std::make_pair(operand_name,change));
    return this->size();
  }
  const unsigned& id(void) const { return id_; }
  const std::string& name(void) const { return name_; }
  unsigned num_operand(void) const { return size(); }
  std::pair<std::string,QN::value_type> action(const unsigned& i) const 
    { return at(i); }
  //const std::string& matrix_element(void) const { return matrixelem_; }
  double matrix_element(const std::string& qn_name, const QN::value_type& qn_state);
private:
  static unsigned count_;
  unsigned id_;
  std::string name_;
  std::string matrixelem_;
};



} // end namespace model

#endif

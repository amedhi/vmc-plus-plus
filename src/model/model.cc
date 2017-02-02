/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-01 22:25:28
*----------------------------------------------------------------------------*/
#include "model.h"

namespace model {

/*unsigned Hamiltonian::add_sitebasis(SiteBasis& sitebasis)
{
  // it's an error if any 'sitebasis' was already added
  if (basis_.size()>0) 
    throw std::logic_error("Hamiltonian::add_sitebasis: 'sitebasis' already exists, overwrite not allowed.");
  // the 'sitebasis' is implicitly defined for all site types
  for (const auto& elem : sitetypes_map_) {
    unsigned mapped_type = elem.second;
    basis_.add_sitebasis(mapped_type,sitebasis); 
  }
  return basis_.size();
}

unsigned Hamiltonian::add_sitebasis(const unsigned& type, SiteBasis& sitebasis)
{
  // add 'sitebasis' of the given 'site type'
  auto it=sitetypes_map_.find(type);
  if (it==sitetypes_map_.end()) 
    throw std::range_error("Hamiltonian::add_sitebasis: specified 'site type' not found");
  unsigned mapped_type = it->second;
  // it's an error if any 'sitebasis' was already added
  if (!basis_.add_sitebasis(mapped_type,sitebasis)) 
    throw std::logic_error("Hamiltonian::add_sitebasis: 'sitebasis' already exists, overwrite not allowed.");
  return basis_.size();
}
*/


unsigned Hamiltonian::add_siteterm(const std::string& name, const CouplingConstant& cc,
  const qn_op& op)
{
  // remap site type values in 'cc'
  CouplingConstant cc_remapped = cc;
  cc_remapped.clear_map();
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    // the 'cc' is implicitly defined for all types
    std::string cc_expr = cc.begin()->second;
    for (const auto& m : sitetypes_map_) {
      cc_remapped.insert({m.second, cc_expr});
    }
  }
  else {
    for (auto it=cc.begin(); it!=cc.end(); ++it) {
      unsigned sitetype = it->first;
      auto it2=sitetypes_map_.find(sitetype);
      if (it2!=sitetypes_map_.end()) {
        unsigned mapped_type = it2->second;
        //std::cout << mapped_type << it->second << "\n";
        cc_remapped.insert({mapped_type, it->second});
      }
      else throw std::range_error("Hamiltonian::add_siteterm: non-existent 'site type' specified");
    }
  }
  unsigned num_sitetypes = sitetypes_map_.size();
  this->std::vector<SiteTerm>::push_back(SiteTerm(name, cc_remapped, op, num_sitetypes));
  return this->std::vector<SiteTerm>::size();
}

unsigned Hamiltonian::add_bondterm(const std::string& name, const CouplingConstant& cc, const qn_op& op)
{
  // remap bond type values in 'cc'
  CouplingConstant cc_remapped = cc;
  cc_remapped.clear_map();
  if (cc.size()==1 && cc.begin()->first==CouplingConstant::global_type) {
    // the 'cc' is implicitly defined for all types
    std::string cc_expr = cc.begin()->second;
    for (const auto& m : bondtypes_map_) {
      //std::cout << "m= " << m.second << " cc_expr = " << cc_expr << "\n";
      cc_remapped.insert({m.second, cc_expr});
    }
  }
  else {
    for (auto it=cc.begin(); it!=cc.end(); ++it) {
      unsigned bondtype = it->first;
      auto it2=bondtypes_map_.find(bondtype);
      if (it2!=bondtypes_map_.end()) {
        unsigned mapped_type = it2->second;
        cc_remapped.insert({mapped_type, it->second});
      }
      else throw std::range_error("Hamiltonian::add_bondterm: non-existent 'site type' specified");
    }
  }
  unsigned num_bondtypes = bondtypes_map_.size();
  std::vector<BondTerm>::push_back(BondTerm(name, cc_remapped, op, num_bondtypes));
  return std::vector<BondTerm>::size();
}

int Hamiltonian::finalize(const lattice::Lattice& L)
{
  // check if 'sitebasis' for all 'site types' are defined
  /*for (const auto& elem : sitetypes_map_) {
    unsigned site_type = elem.second;
    if (basis_.find(site_type) == basis_.end()) 
      throw std::range_error("modellibrary:: 'sitebasis' not defined for all 'site type'-s");
  }*/

  // finalize the site terms
  for (auto it=std::vector<SiteTerm>::begin(); it!=std::vector<SiteTerm>::end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  has_siteterm_ = (std::vector<SiteTerm>::size()>0);
  st_begin_ = std::vector<SiteTerm>::cbegin();
  st_end_ = std::vector<SiteTerm>::cend();
  // finalize the bond terms
  for (auto it=std::vector<BondTerm>::begin(); it!=std::vector<BondTerm>::end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  has_bondterm_ = (std::vector<BondTerm>::size()>0);
  bt_begin_ = std::vector<BondTerm>::cbegin();
  bt_end_ = std::vector<BondTerm>::cend();


  set_info_string(L);
  return 0;
}

void Hamiltonian::change_parameter_value(const std::string& pname, const double& pval) 
{
  auto it = parms_.find(pname);
  if (it != parms_.end()) it->second = pval;
}

void Hamiltonian::update_parameters(const input::Parameters& inputs)
{
  // update the parameter values
  for (auto& p : parms_) p.second = inputs.set_value(p.first, p.second);
  // update the model term couping constants
  for (auto it=std::vector<SiteTerm>::begin(); it!=std::vector<SiteTerm>::end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
  for (auto it=std::vector<BondTerm>::begin(); it!=std::vector<BondTerm>::end(); ++it) {
    it->eval_coupling_constant(constants_, parms_); 
  }
}

double Hamiltonian::get_parameter_value(const std::string& pname) const
{
  auto it = parms_.find(pname);
  if (it != parms_.end()) return it->second;
  else return 0.0;
}

void Hamiltonian::get_term_names(std::vector<std::string>& term_names) const
{
  term_names.clear();
  for (auto it=std::vector<BondTerm>::cbegin(); it!= std::vector<BondTerm>::cend(); ++it) 
    term_names.push_back(it->name());
  for (auto it=std::vector<SiteTerm>::cbegin(); it!= std::vector<SiteTerm>::cend(); ++it) 
    term_names.push_back(it->name());
}

void Hamiltonian::set_info_string(const lattice::Lattice& L) 
{
  info_str_.clear();
  info_str_ << "# Lattice: " << L.name() << " (";
  info_str_ << "Size="<<L.size1()<<"x"<<L.size2()<<"x"<< L.size3()<<", ";
  info_str_ << "Sites/unitcell="<<L.num_basis_sites()<<", ";
  info_str_ << "Boundary="<<static_cast<int>(L.bc1_periodicity()) << "-"; 
  info_str_ << static_cast<int>(L.bc2_periodicity()) << "-";
  info_str_ << static_cast<int>(L.bc3_periodicity()) << ")\n";
  info_str_ << "# No of sites = " << L.num_sites() << "\n";
  info_str_ << "# Hamiltonian: " << model_name << "\n";
  info_str_.precision(6);
  info_str_.setf(std::ios_base::fixed);
  for (const auto& p : parms_) 
    info_str_ << "# " << p.first << " = " << p.second << "\n";
}





} // end namespace model

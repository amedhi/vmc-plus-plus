/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace mc {

Observables::Observables() 
  : energy_("Energy")
  , energy_sq_("Energy^2")
  , magn_("Magnetization")
  , magn_sq_("Magnetization^2")
  , potts_magn_("PottsMagn")
  , potts_magn_sq_("PottsMagn^2")
  , strain_("Strain")
  , strain_sq_("Strain^2")
  , energy_terms_("EnergyTerms")
  , energy_terms_sq_("EnergyTerms^2")
{
  push_back(energy_);
  push_back(energy_sq_);
  push_back(magn_);
  push_back(magn_sq_);
  push_back(potts_magn_);
  push_back(potts_magn_sq_);
  push_back(strain_);
  push_back(strain_sq_);
  push_back(energy_terms_);
  push_back(energy_terms_sq_);
}

void Observables::init(input::Parameters& parms, const model::Model& model,
  void (&print_copyright)(std::ostream& os))
{
  // energy_terms details 
  std::vector<std::string> elem_names;
  model.get_term_names(elem_names);
  energy_terms_.set_elements(elem_names.size(), elem_names);
  // energy_terms_sq details 
  for (auto& name : elem_names) name += "2";
  energy_terms_sq_.set_elements(elem_names.size(), elem_names);

  // actual init of observables
  for (auto& obs : *this) {
    obs.get().init(parms);
    if (obs.get().is_on() && obs.get().is_open() && obs.get().replace_mode()) {
      print_copyright(obs.get().fs_stream());
      model.print_info(obs.get().fs_stream());
      //obs.get().print_heading("T");
    }
  }

  // basis for building operator matrices
  // basis_ = model.basis();
  // build observable operators
  /*if (magn_ || magn_sq_) {
    magn_op_.init(model.basis(), "S(i)");
  }
  if (potts_magn_ || potts_magn_sq_) {
    potts_magn_op_.init(model.basis(), "S(i)");
  }
  if (strain_ || strain_sq_) {
    strain_op_.init(model.basis(), "sigma(i)");
  }*/
}

void Observables::as_function_of(const std::map<std::string, double> xparms)
{
  num_xparms_ = xparms.size();
  for (auto& obs : *this) {
    if (obs.get().is_on() && obs.get().is_open() && obs.get().replace_mode()) {
      obs.get().print_heading(xparms);
    }
  }
} 

void Observables::as_function_of(const std::string& xparm_name)
{
  single_xparm_.clear(); 
  single_xparm_[xparm_name] = 0.0;
  as_function_of(single_xparm_);
}

void Observables::print(const std::map<std::string, double> xparms) 
{
  if (xparms.size() != num_xparms_) 
    throw std::range_error("Observables::print: 'x-parameters' size mismatch");
  for (auto& obs : *this) obs.get().print_result(xparms); 
}

void Observables::print(const double& xparm_val) 
{
  if (num_xparms_ != 1) 
    throw std::range_error("Observables::print: no 'x-parameter' was set initially");
  single_xparm_.begin()->second = xparm_val;
  for (auto& obs : *this) obs.get().print_result(single_xparm_); 
}

int ScalarObservable::print_heading(const std::map<std::string, double> xparms) 
{
  if (!is_on()) return 1;
  if (!fs_) throw std::runtime_error("ScalarObservable::print_heading: file not open");
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << std::left;
  fs_ << "# ";
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xparms) 
    fs_ << std::setw(14)<<p.first.substr(0,14);
  std::string short_name = name();
  if (short_name.size() >= 16) short_name.erase(14);
  fs_ << std::setw(14)<<short_name<<std::setw(11)<<"err";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  return 0;
}

int ScalarObservable::print_result(const std::map<std::string, double> xparms) 
{
  if (!is_on()) return 1;
  if (!fs_) throw std::runtime_error("ScalarObservable::print: file not open");
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xparms) 
    fs_ << std::setw(13) << p.second;
  fs_ << data_.result_str() << data_.conv_str();
  fs_ << std::endl;
  return 0;
} 

int VectorObservable::print_heading(const std::map<std::string, double> xparms) 
{
  if (!is_on()) return 1;
  if (!fs_) throw std::runtime_error("VectorObservable::print_heading: file not open");
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xparms) fs_ << std::setw(14)<<p.first.substr(0,14);
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  fs_ << std::setw(9)<<"samples";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  return 0;
}

int VectorObservable::print_result(const std::map<std::string, double> xparms) 
{
  if (!is_on()) return 1;
  if (!fs_) throw std::runtime_error("VectorObservable::print: file not open");
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xparms) 
    fs_ << std::setw(13) << p.second;
  for (unsigned i=0; i<size_; ++i) 
    fs_ << data_.result_str(i); 
  fs_ << data_.conv_str(0).substr(0,10); 
  fs_ << std::endl;
  return 0;
} 

void ObservableBase::init(input::Parameters& parms)
{
  int no_warn;
  onoff_ = parms.set_value(name(), false, no_warn);
  if (onoff_) {
    std::string fname = name_;
    boost::to_lower(fname);
    auto pos = fname.find('^');
    if (pos != std::string::npos) fname.erase(pos,1);
    fname = "res_"+fname+".txt";
    std::string mode = parms.set_value("mode", "NEW");
    boost::to_upper(mode);
    if (mode=="APPEND") replace_mode_ = false;
    else replace_mode_ = true;
    // open file
    if (replace_mode_) fs_.open(fname);
    else fs_.open(fname, std::ios::app);
  }
} 

void ScalarObservable::init(input::Parameters& parms)
{
  ObservableBase::init(parms);
  if (is_on()) data_.init(name());
} 

void VectorObservable::init(input::Parameters& parms)
{
  if (size_==0) 
    throw std::logic_error("VectorObservable::init: can't initialize zero size observable");
  ObservableBase::init(parms);
  if (is_on()) data_.init(name(),size_);
} 


} // end namespace mc


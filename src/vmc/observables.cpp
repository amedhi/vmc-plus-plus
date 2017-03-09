/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace vmc {

ObservableSet::ObservableSet() 
  : energy_("Energy")
  , total_energy_("TotalEnergy")
  , energy_grad_("EnergyGradient")
  , energy_grad2_("EnergyGradient2")
  , sr_coeffs_("SR_Coefficients")
{
  push_back(energy_);
  push_back(total_energy_);
  push_back(energy_grad_);
  push_back(energy_grad2_);
  push_back(sr_coeffs_);
}

void ObservableSet::init(const input::Parameters& inputs, 
    void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model,
    const std::vector<std::string>& varp_names)
{
  // file open mode
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  // check which observables to calculate
  for (auto& obs : *this) obs.get().check_on(inputs, replace_mode_);
  // heading message
  print_copyright(headstream_);
  model.print_info(headstream_);
  num_xvars_ = 0; 

  // Energy (always set up)
  std::vector<std::string> elem_names;
  model.get_term_names(elem_names);
  energy_.set_elements(elem_names);
  energy_.set_have_total();

  unsigned num_varp = varp_names.size();
  // Energy Gradient
  if (energy_grad_) {
    energy_grad_.set_elements(varp_names);
    energy_grad2_.set_elements(2*num_varp);
  }
  if (energy_||total_energy_||energy_grad_) need_energy_ = true;
}

void ObservableSet::as_functions_of(const std::vector<std::string>& xvars)
{
  xvars_ = xvars;
  num_xvars_ = xvars_.size();
}

void ObservableSet::as_functions_of(const std::string& xvar)
{
  xvars_ = {xvar};
  num_xvars_ = 1;
}

void ObservableSet::switch_off(void) {
  for (auto& obs : *this) obs.get().switch_off();
}

void ObservableSet::print_heading(void)
{
  for (auto& obs : *this) {
    if (obs.get().is_on() && obs.get().replace_mode()) 
      obs.get().print_heading(headstream_,xvars_);
  }
}

void ObservableSet::print_results(const std::vector<double>& xvals) 
{
  if (num_xvars_ == xvals.size()) 
    throw std::invalid_argument("Observables::print_result: 'x-vars' size mismatch");
  for (auto& obs : *this) {
    if (obs.get().is_on()) {
      obs.get().print_heading(headstream_,xvars_);
      obs.get().print_result(xvals);
    }
  }
}

void ObservableSet::print_results(const double& xval) 
{
  if (num_xvars_ != 1) 
    throw std::invalid_argument("Observables::print_result: 'x-vars' size mismatch");
  std::vector<double> xvals{xval};
  for (auto& obs : *this) {
    if (obs.get().is_on()) {
      obs.get().print_heading(headstream_,xvars_);
      obs.get().print_result(xvals);
    }
  }
}

//--------------------Observable class-----------------------
Observable::Observable() 
  : MC_Data() 
{
  num_dataset_=0;
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
}

Observable::Observable(const std::string& name, const unsigned& size) 
{
  this->init(name,size);
} 

void Observable::init(const std::string& name, const unsigned& size) 
{
  MC_Data::init(name, size);
  name_ = name;
  num_dataset_=0;
  avg_mcdata_.init(name, size);
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
  elem_names_ = {name};
  // file name
  fname_ = name_;
  boost::to_lower(fname_);
  auto pos = fname_.find('^');
  if (pos != std::string::npos) fname_.erase(pos,1);
  fname_ = "res_"+fname_+".txt";
}

void Observable::init(const std::string& name, const bool& replace_mode, 
  const unsigned& size) 
{
  init(name, size);
  is_on_ = true;
  replace_mode_ = replace_mode;
}

void Observable::check_on(const input::Parameters& inputs, const bool& replace_mode) 
{
  int no_warn;
  is_on_ = inputs.set_value(name(), false, no_warn);
  replace_mode_ = replace_mode;
}

void Observable::print_heading(const std::stringstream& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!fs_.is_open()) open_file();
  fs_ << header.rdbuf();
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  // total value
  if (MC_Data::size()>1 && have_total_)
    fs_ << std::setw(14)<<"Total"<<std::setw(11)<<"err";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  heading_printed_ = true;
  close_file();
}

void Observable::print_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return;
  if (!fs_.is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  // total value
  //if (MC_Data::size()>1)
  if (MC_Data::size()>1 && have_total_)
    fs_ << MC_Data::result_str(-1); 
  for (unsigned i=0; i<MC_Data::size(); ++i) 
    fs_ << MC_Data::result_str(i); 
  fs_ << MC_Data::conv_str(0); //.substr(0,10); 
  fs_ << std::endl;
  close_file();
} 

void Observable::open_file(void) 
{
  if (fs_.is_open()) return;
  if (replace_mode_) {
    fs_.open(fname_);
    replace_mode_ = false;
  }
  else fs_.open(fname_, std::ios::app);
  if (!fs_.is_open()) 
    throw std::runtime_error("Observable::open_file: file open failed");
}

void Observable::close_file(void) 
{
  fs_.close();
}

void Observable::resize(const unsigned& size) 
{
  MC_Data::resize(size);
  avg_mcdata_.resize(size);
  elem_names_.resize(size);
}

void Observable::set_elements(const std::vector<std::string>& elem_names)
{
  MC_Data::resize(elem_names.size());
  avg_mcdata_.resize(elem_names.size());
  elem_names_ = elem_names;
}

void Observable::set_elements(const unsigned& size)
{
  MC_Data::resize(size);
  avg_mcdata_.resize(size);
  elem_names_.resize(size); 
  for (unsigned i=0; i<size; ++i) 
    elem_names_[i] = "var" + std::to_string(i);
}

void Observable::save_result(void)
{
  avg_mcdata_ << MC_Data::mean_data();
  avg_stddev_ += MC_Data::stddev_data();
  avg_tau_ += MC_Data::tau();
  ++num_dataset_;
}


} // end namespace vmc


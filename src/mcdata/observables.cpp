/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace obs {

ObservableSet::ObservableSet() 
  : energy_("Energy")
  , total_energy_("TotalEnergy")
  , energy_grad_("EnergyGradient")
{
  push_back(energy_);
  push_back(total_energy_);
  push_back(energy_grad_);
}

void ObservableSet::init(const input::Parameters& inputs) 
{
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  for (auto& obs : *this) obs.get().check_on(inputs, replace_mode_);
}

void ObservableSet::switch_off(void) {
  for (auto& obs : *this) obs.get().switch_off();
}

void ObservableSet::file_init(const std::vector<std::string>& xpnames, 
   void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model)
{
  num_xparms_ = xpnames.size();
  for (auto& obs : *this) {
    if (obs.get().is_on() && obs.get().replace_mode()) {
      obs.get().open_file();
      print_copyright(obs.get().fs_stream());
      model.print_info(obs.get().fs_stream());
      obs.get().print_heading(xpnames);
    }
  }
}

void ObservableSet::print_heading(const std::vector<std::string>& xpnames, 
   void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model)
{
  file_init(xpnames,print_copyright,model);
}

void ObservableSet::print_heading(const std::string& xpname, 
   void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model)
{
  std::vector<std::string> xpnames({xpname});
  file_init(xpnames,print_copyright,model);
}

void ObservableSet::print_heading(void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model)
{
  std::vector<std::string> xpnames;
  file_init(xpnames,print_copyright,model);
}

void ObservableSet::print(const std::vector<double> xpvals) 
{
  if (xpvals.size() != num_xparms_) 
    throw std::invalid_argument("Observables::print: 'x-parameters' size mismatch");
  for (auto& obs : *this) obs.get().print_result(xpvals); 
}

void ObservableSet::print(const double& xparm_val) 
{
  if (num_xparms_ != 1) 
    throw std::invalid_argument("Observables::print: no 'x-parameter' was set initially");
  std::vector<double> xpvals({xparm_val});
  for (auto& obs : *this) obs.get().print_result(xpvals); 
}

void ObservableSet::print(void) 
{
  if (num_xparms_ != 0) 
    throw std::invalid_argument("Observables::print: 'x-parameter' value not given");
  std::vector<double> xpvals;
  for (auto& obs : *this) obs.get().print_result(xpvals); 
}

//--------------------Observable class-----------------------
Observable::Observable() 
  : mcdata() 
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
  mcdata::init(name, size);
  name_ = name;
  num_dataset_=0;
  avg_mcdata_.init(name, size);
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
  elem_names_ = {name};
}

void Observable::resize(const unsigned& size) 
{
  mcdata::resize(size);
  avg_mcdata_.resize(size);
  elem_names_.resize(size);
}

void Observable::set_elements(const std::vector<std::string>& elem_names)
{
  mcdata::resize(elem_names.size());
  avg_mcdata_.resize(elem_names.size());
  elem_names_ = elem_names;
}

void Observable::set_elements(const unsigned& size)
{
  mcdata::resize(size);
  avg_mcdata_.resize(size);
  elem_names_.resize(size); 
  for (unsigned i=0; i<size; ++i) 
    elem_names_[i] = "var" + std::to_string(i);
}

void Observable::save_result(void)
{
  avg_mcdata_ << mcdata::mean_data();
  avg_stddev_ += mcdata::stddev_data();
  avg_tau_ += mcdata::tau();
  ++num_dataset_;
}

void Observable::check_on(const input::Parameters& inputs, const bool& replace_mode) 
{
  int no_warn;
  is_on_ = inputs.set_value(name(), false, no_warn);
  replace_mode_ = replace_mode;
}

void Observable::open_file(void) 
{
  if (fs_.is_open()) return;
  if (is_on_) {
    std::string fname = name_;
    boost::to_lower(fname);
    auto pos = fname.find('^');
    if (pos != std::string::npos) fname.erase(pos,1);
    fname = "res_"+fname+".txt";
    // open file
    if (replace_mode_) fs_.open(fname);
    else fs_.open(fname, std::ios::app);
  }
}

int Observable::print_heading(const std::vector<std::string>& xpnames) 
{
  if (!is_on()) return 1;
  open_file();
  if (!fs_.is_open()) throw std::runtime_error("Observable::print_heading: file not open");
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xpnames) fs_ << std::setw(14)<<p.substr(0,14);
  // total value
  if (mcdata::size()>1)
    fs_ << std::setw(14)<<"Total"<<std::setw(11)<<"err";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  return 0;
}

int Observable::print_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return 1;
  open_file();
  if (!fs_.is_open()) throw std::runtime_error("Observable::print_result: file not open");
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  // total value
  if (mcdata::size()>1)
    fs_ << mcdata::result_str(-1); 
  for (unsigned i=0; i<mcdata::size(); ++i) 
    fs_ << mcdata::result_str(i); 
  fs_ << mcdata::conv_str(0); //.substr(0,10); 
  fs_ << std::endl;
  return 0;
} 



} // end namespace mc


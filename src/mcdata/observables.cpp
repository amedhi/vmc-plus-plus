/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace mc {

Observables::Observables(const input::Parameters& inputs) 
  : energy_("Energy")
  , total_energy_("TotalEnergy")
  , energy_grad_("EnergyGradient")
  //, magn_("Magnetization")
  //, magn_sq_("Magnetization^2")
{
  push_back(energy_);
  push_back(total_energy_);
  push_back(energy_grad_);
  //push_back(energy_sq_);
  //push_back(magn_);
  //push_back(magn_sq_);
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  for (auto& obs : *this) obs.get().check_on(inputs, replace_mode_);
}

void Observables::switch_off(void) {
  for (auto& obs : *this) obs.get().switch_off();
}

void Observables::init(const input::Parameters& inputs) {
  bool replace_mode;
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode = false;
  else replace_mode = true;
  for (auto& obs : *this) obs.get().check_on(inputs, replace_mode);
}

void Observables::print_heading(const std::vector<std::string>& xpnames, 
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

void Observables::print_heading(const std::string& xpname, 
   void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model)
{
  num_xparms_ = 1;
  std::vector<std::string> xpnames({xpname});
  for (auto& obs : *this) {
    if (obs.get().is_on() && obs.get().replace_mode()) {
      obs.get().open_file();
      print_copyright(obs.get().fs_stream());
      model.print_info(obs.get().fs_stream());
      obs.get().print_heading(xpnames);
    }
  }
}

void Observables::print_heading(void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model)
{
  num_xparms_ = 0;
  std::vector<std::string> xpnames;
  for (auto& obs : *this) {
    if (obs.get().is_on() && obs.get().replace_mode()) {
      obs.get().open_file();
      print_copyright(obs.get().fs_stream());
      model.print_info(obs.get().fs_stream());
      obs.get().print_heading(xpnames);
    }
  }
}

/*
void Observables::as_function_of(const std::vector<std::string>& xpnames)
{
  num_xparms_ = xpnames.size();
  for (auto& obs : *this) {
    if (obs.get().is_on() && obs.get().is_open() && obs.get().replace_mode()) {
      obs.get().print_heading(xpnames);
    }
  }
} 
void Observables::as_function_of(const std::string& xpname)
{
  std::vector<std::string> xpnames({xpname});
  as_function_of(xpnames);
}
*/

void Observables::print(const std::vector<double> xpvals) 
{
  if (xpvals.size() != num_xparms_) 
    throw std::invalid_argument("Observables::print: 'x-parameters' size mismatch");
  for (auto& obs : *this) obs.get().print_result(xpvals); 
}

void Observables::print(const double& xparm_val) 
{
  if (num_xparms_ != 1) 
    throw std::invalid_argument("Observables::print: no 'x-parameter' was set initially");
  std::vector<double> xpvals({xparm_val});
  for (auto& obs : *this) obs.get().print_result(xpvals); 
}

void Observables::print(void) 
{
  if (num_xparms_ != 0) 
    throw std::invalid_argument("Observables::print: 'x-parameter' value not given");
  std::vector<double> xpvals;
  for (auto& obs : *this) obs.get().print_result(xpvals); 
}

int ScalarObservable::print_heading(const std::vector<std::string>& xpnames) 
{
  if (!is_on()) return 1;
  open_file();
  if (!fs_.is_open()) throw std::runtime_error("* ScalarObservable::print_heading: file not open");
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << std::left;
  fs_ << "# ";
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xpnames) 
    fs_ << std::setw(14)<<p.substr(0,14);
  std::string short_name = name();
  if (short_name.size() >= 16) short_name.erase(14);
  fs_ << std::setw(14)<<short_name<<std::setw(11)<<"err";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  return 0;
}

int ScalarObservable::print_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return 1;
  open_file();
  if (!fs_.is_open()) throw std::runtime_error("ScalarObservable::print: file not open");
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  fs_ << data_.result_str() << data_.conv_str();
  fs_ << std::endl;
  return 0;
} 

int VectorObservable::print_heading(const std::vector<std::string>& xpnames) 
{
  if (!is_on()) return 1;
  open_file();
  if (!fs_.is_open()) throw std::runtime_error("VectorObservable::print_heading: file not open");
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xpnames) fs_ << std::setw(14)<<p.substr(0,14);
  // total value
  fs_ << std::setw(14)<<"Total"<<std::setw(11)<<"err";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  fs_ << std::setw(9)<<"samples";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  return 0;
}

int VectorObservable::print_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return 1;
  open_file();
  if (!fs_.is_open()) throw std::runtime_error("* VectorObservable::print: file not open");
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  // total value
  fs_ << data_.result_str(-1); 
  for (unsigned i=0; i<size_; ++i) 
    fs_ << data_.result_str(i); 
  fs_ << data_.conv_str(0).substr(0,10); 
  fs_ << std::endl;
  return 0;
} 

void ObservableBase::check_on(const input::Parameters& inputs, const bool& replace_mode)
{
  int no_warn;
  is_on_ = inputs.set_value(name(), false, no_warn);
  replace_mode_ = replace_mode;
}

void ScalarObservable::check_on(const input::Parameters& inputs, const bool& replace_mode)
{
  ObservableBase::check_on(inputs,replace_mode);
  if (is_on()) data_.init(name());
}

void VectorObservable::check_on(const input::Parameters& inputs, const bool& replace_mode)
{
  ObservableBase::check_on(inputs,replace_mode);
}

void ScalarObservable::switch_on(void) 
{
  ObservableBase::switch_on();
  data_.init(name());
}

void VectorObservable::switch_on(void)
{
  ObservableBase::switch_on(); 
  size_ = 0;
}
/*void VectorObservable::switch_on(const unsigned& size)
{
  size_=size;   
  if (size_==0) 
    throw std::logic_error("VectorObservable::switch_on: can't switch on zero size observable");
  ObservableBase::switch_on(); 
  elem_names_=std::vector<std::string>(size_); 
  data_.init(name(),size_);
}*/

void ObservableBase::open_file(void) 
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


void VectorObservable::set_elements(const std::vector<std::string>& elem_names)
{
  elem_names_=elem_names; 
  size_=elem_names_.size();   
  if (size_==0) 
    throw std::logic_error("VectorObservable::set_elements: can't initialize zero size observable");
  data_.init(name(),size_);
}


} // end namespace mc


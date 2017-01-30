/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef MC_OBSERVABLES_H
#define MC_OBSERVABLES_H

#include <string>
#include <vector>
#include <map>
#include <functional>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include <Eigen/Core>
#include "../scheduler/inputparams.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "mcdata.h"
//#include "observable_operator.h"

namespace mc {

using RealObservableData = mcdata<double>;
using VectorObservableData = mcdata<VectorData>;

class ObservableBase 
{
public:
  ObservableBase(const std::string& name) : name_{name} {}
  ~ObservableBase() { fs_.close(); }
  virtual void init(input::Parameters& parms); 
  virtual void reset(void) {}
  virtual int print_heading(const std::map<std::string, double> xparms) { return -1; }
  virtual int print_result(const std::map<std::string, double> xparms) { return -1; }
  inline const bool& is_on(void) const { return onoff_; }
  bool& state(void) { return onoff_; }
  void open_file(const std::string& fname) { if (onoff_) fs_.open(fname); }
  bool is_open(void) const { if (fs_) return true; else return false; }
  const bool& replace_mode(void) const { return replace_mode_; }
  std::ofstream& fs_stream(void) { return fs_; }  
  const std::string& name(void) const { return name_; }
protected:
  std::ofstream fs_;
private:
  bool onoff_{false};
  bool replace_mode_{true};
  std::string name_{""};
};

class ScalarObservable : public ObservableBase 
{
public:
  ScalarObservable(const std::string& name) : ObservableBase(name) {}
  ~ScalarObservable() {}
  void init(input::Parameters& parms) override; 
  inline void operator<<(const double& data) { data_ << data; }
  int print_heading(const std::map<std::string, double> xparms) override;
  int print_result(const std::map<std::string, double> xparms) override;
  void reset(void) override { if (is_on()) data_.clear(); }
  operator int(void) const { return is_on(); }
private:
  RealObservableData data_;
};

class VectorObservable : public ObservableBase {
public:
  VectorObservable(const std::string& name) : ObservableBase(name) {}
  //VectorObservable(const std::string& name, const std::vector<std::string>& elem_names, 
  //  const unsigned& size) { init(name, elem_names, size); } 
  ~VectorObservable() {}
  void set_elements(const unsigned& size, const std::vector<std::string>& elem_names)
    { size_=size; elem_names_=elem_names; }
  void init(input::Parameters& parms) override; 
  inline void operator<<(const VectorData& data) { data_ << data; }
  int print_heading(const std::map<std::string, double> xparms) override;
  int print_result(const std::map<std::string, double> xparms) override;
  void reset(void) override { if (is_on()) data_.clear(); }
  operator int(void) const { return is_on(); }
  const unsigned& size(void) const { return size_; }
private:
  unsigned size_{0};
  VectorObservableData data_;
  std::vector<std::string> elem_names_;
};

class Observables : public std::vector<std::reference_wrapper<ObservableBase> >
{
public:
  Observables();
  ~Observables() {}
  void init(input::Parameters& parms, const model::Model& model, 
    void (&print_copyright)(std::ostream& os));
  inline ScalarObservable& energy(void) { return energy_; }
  inline ScalarObservable& energy_sq(void) { return energy_sq_; }
  inline ScalarObservable& magn(void) { return magn_; }
  inline ScalarObservable& magn_sq(void) { return magn_sq_; }
  inline ScalarObservable& potts_magn(void) { return potts_magn_; }
  inline ScalarObservable& potts_magn_sq(void) { return potts_magn_sq_; }
  inline ScalarObservable& strain(void) { return strain_; }
  inline ScalarObservable& strain_sq(void) { return strain_sq_; }
  inline VectorObservable& energy_terms(void) { return energy_terms_; }
  inline VectorObservable& energy_terms_sq(void) { return energy_terms_sq_; }
  void reset(void) { for (auto& obs : *this) obs.get().reset(); }
  //inline SiteObsOperator& magn_op(void) { return magn_op_; } 
  //inline SiteObsOperator& potts_magn_op(void) { return potts_magn_op_; } 
  //inline SiteObsOperator& strain_op(void) { return strain_op_; } 
  //void set_magn_op(const std::string& op, const std::string& site="i"); 
  //void set_strain_op(const std::string& op, const std::string& site="i"); 
  //const bool& need_energy(void) { return energy_.on(); }
  void as_function_of(const std::map<std::string, double> xparms);
  void as_function_of(const std::string& xparm_name);
  void print(const std::map<std::string, double> xparms); 
  void print(const double& xparm_val);
private:
  ScalarObservable energy_;
  ScalarObservable energy_sq_;
  ScalarObservable magn_;
  ScalarObservable magn_sq_;
  ScalarObservable potts_magn_;
  ScalarObservable potts_magn_sq_;
  ScalarObservable strain_;
  ScalarObservable strain_sq_;
  VectorObservable energy_terms_;
  VectorObservable energy_terms_sq_;
  //SiteObsOperator magn_op_;
  //SiteObsOperator potts_magn_op_;
  //SiteObsOperator strain_op_;
  //model::BasisDescriptor basis_;
  unsigned num_xparms_{0};
  std::map<std::string, double> single_xparm_;
};



} // end namespace mc

#endif

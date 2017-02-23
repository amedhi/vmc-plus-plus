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
#include "../scheduler/task.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "./mcdata.h"

namespace mc {

using RealObservableData = mcdata<double>;
using VectorObservableData = mcdata<VectorData>;

class ObservableBase 
{
public:
  ObservableBase(const std::string& name) : name_{name} {}
  ~ObservableBase() { fs_.close(); }
  virtual void check_on(const input::Parameters& parms, const bool& replace_mode); 
  virtual void reset(void) {}
  virtual int print_heading(const std::vector<std::string>& xpnames) { return -1; }
  virtual int print_result(const std::vector<double>& xpvals) { return -1; }
  virtual void switch_on(void) { is_on_=true; }
  void switch_off(void) { is_on_=false; if (fs_.is_open()) fs_.close(); }
  void open_file(void); 
  inline const bool& is_on(void) const { return is_on_; }
  bool is_open(void) const { return fs_.is_open(); }
  const bool& replace_mode(void) const { return replace_mode_; }
  std::ofstream& fs_stream(void) { return fs_; }  
  const std::string& name(void) const { return name_; }
protected:
  std::ofstream fs_;
private:
  bool is_on_{false};
  bool replace_mode_{true};
  std::string name_{""};
};

class ScalarObservable : public ObservableBase 
{
public:
  ScalarObservable(const std::string& name) : ObservableBase(name) {}
  ~ScalarObservable() {}
  virtual void check_on(const input::Parameters& parms, const bool& replace_mode) override; 
  inline void operator<<(const double& data) { data_ << data; }
  void switch_on(void) override; 
  int print_heading(const std::vector<std::string>& xpnames) override;
  int print_result(const std::vector<double>& xpvals) override;
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
  virtual void check_on(const input::Parameters& parms, const bool& replace_mode) override; 
  void switch_on(void) override; 
  //void switch_on(const unsigned& size); 
  void set_elements(const std::vector<std::string>& elem_names);
  inline void operator<<(const VectorData& data) { data_ << data; }
  int print_heading(const std::vector<std::string>& xpnames) override;
  int print_result(const std::vector<double>& xpvals) override;
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
  Observables(const input::Parameters& inputs);
  ~Observables() {}
  void switch_off(void);
  void init(const input::Parameters& inputs);
  void print_heading(void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
  void print_heading(const std::string& xpname, void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
  void print_heading(const std::vector<std::string>& xpnames, void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
  //inline ScalarObservable& magn(void) { return magn_; }
  //inline ScalarObservable& magn_sq(void) { return magn_sq_; }
  inline VectorObservable& energy(void) { return energy_; }
  inline ScalarObservable& total_energy(void) { return total_energy_; }
  inline VectorObservable& energy_grad(void) { return energy_grad_; }
  //inline VectorObservable& energy_sq(void) { return energy_sq_; }
  void reset(void) { for (auto& obs : *this) obs.get().reset(); }
  //inline SiteObsOperator& magn_op(void) { return magn_op_; } 
  //inline SiteObsOperator& potts_magn_op(void) { return potts_magn_op_; } 
  //inline SiteObsOperator& strain_op(void) { return strain_op_; } 
  //void as_function_of(const std::vector<std::string>& xpnames);
  //void as_function_of(const std::string& xpname);
  void print(const std::vector<double> xpvals); 
  void print(const double& xpval);
  void print(void);
private:
  bool replace_mode_{true};
  //ScalarObservable magn_;
  //ScalarObservable magn_sq_;
  VectorObservable energy_;
  ScalarObservable total_energy_;
  VectorObservable energy_grad_;
  //VectorObservable energy_sq_;
  unsigned num_xparms_{0};
  //std::map<std::string, double> single_xparm_;
};


} // end namespace mc

#endif

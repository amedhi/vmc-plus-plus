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

namespace obs {

using vector = mcdata::data_t;
using scalar = mcdata::scalardata_t;

class Observable : public mcdata
{
public:
  using data_t = DataBin::data_t;
  Observable();
  Observable(const std::string& name, const unsigned& size=1);
  ~Observable() {}
  void init(const std::string& name, const unsigned& size=1) override; 
  void resize(const unsigned& size) override;
  void set_option(const bool& file_out=true, const bool& print_total=false);
  void set_elements(const std::vector<std::string>& elem_names);
  void set_elements(const unsigned& size);
  void reset(void) { mcdata::clear(); }
  void check_on(const input::Parameters& inputs, const bool& replace_mode); 
  void switch_on(void) { is_on_=true; }
  void switch_off(void) { is_on_=false; if (fs_.is_open()) fs_.close(); }
  operator int(void) const { return is_on(); }
  const bool& is_on(void) const { return is_on_; }
  void open_file(void); 
  bool is_open(void) const { return fs_.is_open(); }
  std::ofstream& fs_stream(void) { return fs_; }  
  const bool& replace_mode(void) const { return replace_mode_; }
  int print_heading(const std::vector<std::string>& xpnames); 
  int print_result(const std::vector<double>& xpvals); 
  void save_result(void);
private:
  unsigned num_dataset_{0};
  mcdata avg_mcdata_;
  data_t avg_stddev_;
  double avg_tau_;
  // for printing
  std::vector<std::string> elem_names_;
  std::ofstream fs_;
  bool is_on_{false};
  bool replace_mode_{true};
  bool need_file_out_{true};
  bool print_total_{false};
  std::string name_{""};
};

class ObservableSet : private std::vector<std::reference_wrapper<Observable> >
{
public:
  ObservableSet();
  ~ObservableSet() {}
  void init(const input::Parameters& inputs);
  void switch_off(void);
  void print_heading(void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
  void print_heading(const std::string& xpname, void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
  void print_heading(const std::vector<std::string>& xpnames, void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
  void reset(void) { for (auto& obs : *this) obs.get().reset(); }
  inline Observable& energy(void) { return energy_; }
  inline Observable& total_energy(void) { return total_energy_; }
  inline Observable& energy_grad(void) { return energy_grad_; }
  inline Observable& energy_grad2(void) { return energy_grad2_; }
  void print(const std::vector<double> xpvals); 
  void print(const double& xpval);
  void print(void);
private:
  bool replace_mode_{true};
  Observable energy_;
  Observable total_energy_;
  Observable energy_grad_;
  Observable energy_grad2_;
  unsigned num_xparms_{0};
  
  void file_init(const std::vector<std::string>& xpnames, void (&print_copyright)(std::ostream& os), 
    const model::Hamiltonian& model);
};


} // end namespace mc

#endif

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
#include <sstream>
#include <iomanip>
#include "../scheduler/task.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../mcdata/mcdata.h"

namespace vmc {

enum class observable_set {normal, energy, energy_grad, sr_coeffs};

class Observable : public mcdata::MC_Data
{
public:
  using data_t = mcdata::data_t;
  using scalar_t = mcdata::scalardata_t;
  Observable();
  Observable(const std::string& name, const unsigned& size=1);
  ~Observable() {}
  void init(const std::string& name, const unsigned& size=1) override; 
  void init(const std::string& name, const bool& replace_mode, const unsigned& size=1); 
  void resize(const unsigned& size) override;
  void save_result(void);
  void set_have_total(void) { have_total_=true; }
  void set_elements(const std::vector<std::string>& elem_names);
  void set_elements(const unsigned& size);
  void reset(void) { MC_Data::clear(); }
  void check_on(const input::Parameters& inputs, const bool& replace_mode); 
  void switch_on(void) { is_on_=true; }
  void switch_on(const unsigned& size) { is_on_=true; set_elements(size); }
  void switch_off(void) { is_on_=false; if (fs_.is_open()) fs_.close(); }
  operator int(void) const { return is_on(); }
  const bool& is_on(void) const { return is_on_; }
  bool is_open(void) const { return fs_.is_open(); }
  std::ofstream& fstream(void) { open_file(); return fs_; }  
  const bool& replace_mode(void) const { return replace_mode_; }
  void print_heading(const std::stringstream& header, 
    const std::vector<std::string>& xvars);
  void print_result(const std::vector<double>& xvals); 
private:
  unsigned num_dataset_{0};
  MC_Data avg_mcdata_;
  data_t avg_stddev_;
  double avg_tau_;
  // for printing
  std::vector<std::string> elem_names_;
  std::ofstream fs_;
  bool is_on_{false};
  bool replace_mode_{true};
  bool heading_printed_{false};
  bool have_total_{false};
  std::string name_{""};
  std::string fname_{""};

  void open_file(void); 
  void close_file(void); 
};

class SC_Correlation : public Observable
{
public:
  using site_t = lattice::LatticeGraph::site_descriptor;
  using Observable::Observable;
  void setup(const lattice::LatticeGraph& graph);
  const unsigned& num_site_pairs(void) const { return src_pairs_size_; }
  const std::pair<site_t,site_t>& site_pair(const unsigned& i) const 
    { return src_pairs_[i]; }
private:
  int max_dist_{0};
  unsigned num_bond_types_{0};
  unsigned src_pairs_size_{0};
  std::vector<std::pair<site_t,site_t> > src_pairs_;
  std::vector<int> pair_distance_;
  std::vector<Eigen::MatrixXd> bond_pair_corr_;
  std::vector<Eigen::MatrixXi> num_symm_pairs_;
};


class ObservableSet : private std::vector<std::reference_wrapper<Observable> >
{
public:
  ObservableSet();
  ~ObservableSet() {}
  void init(const input::Parameters& inputs, 
    void (&print_copyright)(std::ostream& os), const model::Hamiltonian& model,
    const std::vector<std::string>& varp_names);
  void as_functions_of(const std::vector<std::string>& xvars=std::vector<std::string>());
  void as_functions_of(const std::string& xvar);
  void switch_off(void);
  void reset(void) { for (auto& obs : *this) obs.get().reset(); }
  inline Observable& energy(void) { return energy_; }
  inline Observable& total_energy(void) { return total_energy_; }
  inline Observable& energy_grad(void) { return energy_grad_; }
  inline Observable& energy_grad2(void) { return energy_grad2_; }
  inline Observable& sr_coeffs(void) { return sr_coeffs_; }
  inline Observable& sccf(void) { return sccf_; }
  inline SC_Correlation& sc_corr(void) { return sc_corr_; }
  const bool& need_energy(void) const { return need_energy_; }
  void print_heading(void);
  void print_results(const std::vector<double>& xvals=std::vector<double>()); 
  void print_results(const double& xval); 
private:
  bool replace_mode_{true};
  std::stringstream headstream_;
  std::vector<std::string> xvars_;
  unsigned num_xvars_{0};
  Observable energy_;
  Observable total_energy_;
  Observable energy_grad_;
  Observable energy_grad2_;
  Observable sr_coeffs_;
  Observable sccf_;
  SC_Correlation sc_corr_;
  bool need_energy_{false};
};





} // end namespace vmc

#endif

/*---------------------------------------------------------------------------
* InputParameter: Classes for processing input parameters.
* Copyright (C) 2015-2015 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2015-08-17 13:33:19
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-13 11:46:55
*----------------------------------------------------------------------------*/
// File: inputparams.h 

#ifndef INPUT_PARAMETERS_H
#define INPUT_PARAMETERS_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include "./cmdargs.h"

namespace input {

enum class ptype {boo, num, str, nan};
struct pval {bool is_const; ptype type; bool bool_val; double num_val; std::string str_val;};


class JobInput;  // forward declaration

class Parameters : private std::map<std::string, pval> 
{
public:
  friend class JobInput;
  using super_type = std::map<std::string, pval>;
  Parameters() {};
  ~Parameters() {};
  void clear(void) { super_type::clear(); }
  const bool& have_option_quiet(void) const { return have_option_quiet_; }
  const bool& have_option_test(void) const { return have_option_test_; }
  int set_value(const std::string& pname, const int& defval) const;  
  int set_value(const std::string& pname, const int& defval, int& info) const;
  int set_value(const std::string& pname, const int& defval, int& info, bool& is_const) const;
  double set_value(const std::string& pname, const double& defval) const;  
  double set_value(const std::string& pname, const double& defval, int& info) const;  
  double set_value(const std::string& pname, const double& defval, int& info, bool& is_const) const;  
  std::string set_value(const std::string& pname, const std::string& defval) const;  
  std::string set_value(const std::string& pname, const std::string& defval, int& info) const;  
  bool is_constant(const std::string& pname) const;  
  unsigned task_id(void) const { return current_task_; }
  unsigned task_size(void) const { return total_tasks_; }
  void show(const unsigned&) const;
private:
  bool have_option_quiet_{false};
  bool have_option_test_{false};
  std::map<std::string, pval> params;
  unsigned current_task_{0};
  unsigned total_tasks_{0};
  mutable std::map<std::string, pval>::const_iterator it;
  void warn_not_found(const std::string& pname) const;
  void warn_type_mismatch(const std::string& pname, const std::string& type) const;
};

class JobInput 
{
public:
  JobInput() {n_tasks=n_params=0; valid_=false;} 
  JobInput(const scheduler::CommandArg& cmdarg); 
  ~JobInput() {} 
  //bool read_inputs(const std::string& inputfile);
  bool read_inputs(const scheduler::CommandArg& cmdarg);
  //bool read_params(const std::string& inputfile);
  int init_task_params(void);
  int set_task_params(const unsigned& task_id);
  const Parameters& task_params(const unsigned& task_id) 
    { set_task_params(task_id); return task_params_; }
  bool not_valid(void) const {return !valid_;};
  const unsigned& task_size(void) {return n_tasks;}
  //void get_task_param(const unsigned& task_id); 
  //void init_task_parameters(Parameters& p);
  //void set_task_parameters(Parameters& p, const unsigned& task_id); 
private:
  struct parameter {std::string name; ptype type; unsigned size;};
  unsigned int n_params;
  unsigned int n_tasks;
  bool valid_;
  std::string infile;
  std::vector<parameter> param_list;
  std::map<std::string, std::vector<bool> > boo_params;
  std::map<std::string, std::vector<double> > num_params;
  std::map<std::string, std::vector<std::string> > str_params;
  // task parameters
  Parameters task_params_;

  // bad_input exception
  class bad_input: public std::runtime_error
  {
  public:
    bad_input(const std::string& msg, const int& ln=-1);
    std::string message(void) const; 
  private:
    int lnum;
  };
  unsigned int parse(const std::string& inputfile);
};


} // end namespace input

#endif

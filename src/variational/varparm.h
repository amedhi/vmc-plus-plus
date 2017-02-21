/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-20 10:18:07
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-21 10:08:30
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef VARPARM_H
#define VARPARM_H

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include "../scheduler/task.h"

namespace var {

struct parm_t
{
  std::string name;
  double val;
  double lower;
  double upper;
};

class VariationalParms : public std::map<std::string,unsigned>
{
public:
  VariationalParms(); 
  ~VariationalParms() { --num_sets_; }
  void clear(void);
  int add(const std::string& name, const double& value, 
    const double& lower, const double& upper);
  int append_set(const VariationalParms& other_set);
  void update(const input::Parameters& inputs);
  void update(const std::vector<double>& pvec, const unsigned& begin, const unsigned& end);
  const unsigned& id(void) const { return my_id_; }
  parm_t operator[](const unsigned& i) const;
  const std::vector<double>& values(void) const { return val_; }
  const std::vector<double>& lbounds(void) const { return lower_; }
  const std::vector<double>& ubounds(void) const { return upper_; }
  //int update_set(VariationalParms& other_set);
private:

  using idx_pair = std::pair<unsigned,unsigned>;
  std::vector<double> val_;
  std::vector<double> lower_;
  std::vector<double> upper_;
  std::map<unsigned,idx_pair> set_idx_;
  //std::vector<std::pair<> > set_offset_;
  unsigned num_parms_{0};
  unsigned my_id_{0};
  //std::vector<std::pair<it_begin, it_end> > v;
  static unsigned num_sets_; 

  void append(const std::string& name, const double& value, const double& lower, 
    const double& upper);
};

} // end namespace var

#endif
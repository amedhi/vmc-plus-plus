/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef MCDATA_H
#define MCDATA_H

#include <iostream>
#include <sstream>
#include <iomanip>
#include <string>
#include <vector>
#include <stdexcept>
#include <Eigen/Core>

namespace mc {

typedef typename Eigen::Array<double,Eigen::Dynamic,1> VectorData;

/*----------------------DataBin class------------------*/
template<typename T=double>
class DataBin 
{
public:
  using element_type = T;
  DataBin(const int& nrows=1, const int& ncols=1);
  ~DataBin() {}
  void clear(void);
  bool add_sample(const T&);
  //bool add_carried_sample(const T&);
  bool has_samples(void) const { return (num_samples_ > 0); }
  bool has_carry_over(void) const { return !waiting_sample_exist_; }
  const unsigned& num_samples(void) const { return num_samples_; }
  const element_type& carry(void) const { return carry_; } 
  bool have_new_samples(void) const { return num_samples_!=num_samples_last_; }
  void finalize(void) const;
  const T& mean(void) const { finalize(); return mean_; }
  const T& stddev(void) const { finalize(); return stddev_; }
private:
  const element_type& square_root(element_type&) const;
  unsigned num_samples_{0};
  mutable unsigned num_samples_last_{0};
  element_type ssum_;
  element_type sumsq_;
  element_type carry_;
  element_type waiting_sample_;
  bool waiting_sample_exist_{false};
  mutable element_type mean_;
  mutable element_type stddev_;
  element_type Zero_;
  element_type MinusOne_;
  // for vector data
  int num_rows_{1};
  int num_cols_{1};
};

template<typename T>
inline DataBin<T>::DataBin(const int& nrows, const int& ncols) : num_samples_{0}, num_samples_last_{0}, ssum_{0.0}, 
  sumsq_{0.0}, carry_{0.0}, waiting_sample_{0.0}, waiting_sample_exist_{false}, 
  mean_{0.0}, stddev_{-1.0}, Zero_{0.0}, MinusOne_{-1.0}, num_rows_{1}, num_cols_{1}
{
  num_rows_ = nrows;
  num_cols_ = ncols;
}

template<>
inline DataBin<VectorData>::DataBin(const int& nrows, const int& ncols) 
  : num_samples_{0}, num_samples_last_{0}, waiting_sample_exist_{false}
{
  num_rows_ = nrows;
  num_cols_ = ncols;
  Zero_ = VectorData::Zero(num_rows_,num_cols_);
  MinusOne_ = VectorData::Constant(num_rows_,num_cols_,-1.0);
  ssum_ = Zero_;
  sumsq_ = Zero_;
  carry_ = Zero_;
  waiting_sample_ = Zero_;
  mean_ = Zero_;
  stddev_ = MinusOne_;
}

template<typename T>
void DataBin<T>::clear(void) 
{
  num_samples_ = 0;
  num_samples_last_ = 0;
  ssum_ = Zero_;
  sumsq_ = Zero_;
  carry_ = Zero_;
  waiting_sample_ = Zero_;
  waiting_sample_exist_ = false;
  mean_ = Zero_;
  stddev_ = MinusOne_;
}

template<typename T>
bool DataBin<T>::add_sample(const T& new_sample) {
  num_samples_++;
  ssum_ += new_sample;
  sumsq_ += new_sample * new_sample;
  if (waiting_sample_exist_) {
    carry_ = (waiting_sample_ + new_sample)*0.5;
    waiting_sample_exist_ = false;
  }
  else {
    waiting_sample_ = new_sample;
    waiting_sample_exist_ = true;
  }
  return !waiting_sample_exist_;
}

/*template<typename T>
bool  DataBin<T>::add_carried_sample(const T& carried_sample) {
  if (waiting_sample_exist_) {
    // new sample is 'carry_' for next higher bin level 
    carry_ = (waiting_sample_exist_ + carried_sample) * 0.5;
    ssum_ += carry_;
    sumsq_ += carry_ * carry_;
    num_samples_++;
    waiting_sample_exist_ = false;
  }
  else {
    waiting_sample_ = carried_sample;
    waiting_sample_exist_ = true;
  }
  return waiting_sample_exist_;
}
*/

template<typename T>
inline void DataBin<T>::finalize(void) const
{
  if (num_samples_last_ != num_samples_) {
    num_samples_last_ = num_samples_;
    if (num_samples_ > 1) {
      mean_ = ssum_/num_samples_;
      element_type variance_ = sumsq_/num_samples_ - mean_ * mean_;
      stddev_ = variance_/(num_samples_-1);
      stddev_ = square_root(stddev_);
    }
    else {
      mean_ = ssum_;
      stddev_ = MinusOne_;  // to mean infinity
    }
  }
}

template<typename T>
inline const T& DataBin<T>::square_root(element_type& val) const
{
  val = std::sqrt(val);
  return val;
}

template<>
inline const VectorData& DataBin<VectorData>::square_root(VectorData& val) const
{
  val.sqrt();
  return val;
}

/*----------------------mcdata class------------------*/
template<typename T=double>
class mcdata : public std::vector<DataBin<T> >
{
public:
  using element_type = T;
  mcdata() {}
  mcdata(const std::string& name);
  mcdata(const std::string& name, const int& nrows, const int& ncols=1);
  ~mcdata() {}
  void init(const std::string& name);
  void init(const std::string& name, const int& nrows, const int& ncols=1);
  void clear(void);
  void add_sample(const T&);
  void operator<<(const T&);
  const unsigned& num_samples(void) const { return top_bin->num_samples(); }
  void finalize(void) const;
  const element_type& mean(void) const; 
  const element_type& stddev(void) const; // { return top_bin->stddev(); } 
  //const element_type& tau(void) const { finalize(); return tau_; } 
  const std::string& name(void) const { return name_; }
  std::string result_str(void) const; 
  std::string result_str(const int& n) const; 
  std::string conv_str(const int& n=-1) const; 
  const mcdata<T>& with_statistic(void) const { show_statistic_=true; return *this; }
  void show_statistic(std::ostream& os=std::cout) const;
  template<typename S>
  friend std::ostream& operator<<(std::ostream& os, const mcdata<S>& obs);
private:
  void find_conv_and_tau(void) const;
  void find_conv_and_tau(const unsigned& n) const;
  void check_convergence(const std::vector<double>& xv, const std::vector<double>& yv) const;
  std::string name_;
  static const unsigned max_binlevel_defult_ = 20;
  static const unsigned good_sample_size_ = 30;
  typename std::vector<DataBin<T> >::iterator top_bin;
  typename std::vector<DataBin<T> >::iterator end_bin;
  mutable unsigned dcorr_level_;
  mutable element_type mean_;
  mutable element_type stddev_;
  mutable double tau_;
  mutable bool show_statistic_;
  mutable std::string error_converged_;
  mutable std::string convergence_str_;
};

template<typename T>
mcdata<T>::mcdata(const std::string& name) 
{
  this->init(name);
}

template<typename T>
mcdata<T>::mcdata(const std::string& name, const int& nrows, const int& ncols)
{
  this->init(name, nrows, ncols);
}

template<typename T>
void mcdata<T>::init(const std::string& name) 
{
  for (unsigned i=0; i<max_binlevel_defult_; ++i) this->push_back(DataBin<T>());
  top_bin = this->begin();
  end_bin = this->end();
  name_ = name;
  dcorr_level_ = 0;
  mean_ = top_bin->mean();
  stddev_ = top_bin->stddev();
  tau_ = -1.0;
  show_statistic_ = false;
  error_converged_ = "NOT_CONVD";
  convergence_str_ = "NULL";
}

template<typename T>
void mcdata<T>::init(const std::string& name, const int& nrows, const int& ncols) 
{
  for (unsigned i=0; i<max_binlevel_defult_; ++i) this->push_back(DataBin<T>(nrows, ncols));
  top_bin = this->begin();
  end_bin = this->end();
  name_ = name;
  dcorr_level_ = 0;
  mean_ = top_bin->mean();
  stddev_ = top_bin->stddev();
  tau_ = -1.0;
  show_statistic_ = false;
  error_converged_ = "NOT_CONVD";
  convergence_str_ = "NULL";
}

template<typename T>
void mcdata<T>::clear(void) 
{
  for (auto it=this->begin(); it!=this->end(); ++it) it->clear();
  dcorr_level_ = 0;
  mean_ = top_bin->mean();
  stddev_ = top_bin->stddev();
  tau_ = -1.0;
  show_statistic_ = false;
  error_converged_ = "NOT_CONVD";
  convergence_str_ = "NULL";
}

template<typename T>
void  mcdata<T>::operator<<(const T& sample) {
  add_sample(sample);
}

template<typename T>
void mcdata<T>::add_sample(const T& sample) 
{
  auto this_bin = top_bin;  
  element_type new_sample(sample);
  while (this_bin->add_sample(new_sample)) {
    new_sample = this_bin->carry();
    if (this_bin++ == end_bin) break;
  }
}

template<typename T>
const T& mcdata<T>::mean(void) const 
{ 
  this->finalize(); return mean_;
}

template<typename T>
const T& mcdata<T>::stddev(void) const 
{ 
  this->finalize(); return stddev_;
} 

template<typename T>
void mcdata<T>::finalize(void)  const
{ 
  if (top_bin->have_new_samples()) {
    // mean
    mean_ = top_bin->mean();
    // stddev is the 'bin stddev' of the lowest depth bin with a 'good sample size'
    stddev_ = top_bin->stddev();
    dcorr_level_ = 0;
    for (auto bin=this->begin(); bin!=this->end(); ++bin) {
      if (bin->num_samples() < good_sample_size_) break;
      stddev_ = bin->stddev();
      dcorr_level_++;
    }
  }
}

template<typename T>
inline void mcdata<T>::find_conv_and_tau(void) const 
{ 
  this->finalize();
  std::vector<double> xv, yv;
  unsigned level_n = dcorr_level_;
  if (dcorr_level_ == 2) level_n = dcorr_level_-2;
  else if (dcorr_level_ >= 3) level_n = dcorr_level_-3;
  for (unsigned i=level_n; i<=dcorr_level_; ++i) {
    xv.push_back(static_cast<double>(i));
    T stddev = this->operator[](i).stddev();
    yv.push_back(static_cast<double>(stddev));
  }
  // assuming convergence
  double stddev_0 = static_cast<double>(this->operator[](0).stddev());
  double stddev_d = static_cast<double>(this->operator[](dcorr_level_).stddev());
  if (stddev_0 < 1.0E-12) tau_ = -1.0;
  else {
    double r = stddev_d/stddev_0;
    tau_ = 0.5*( r*r - 1.0);
    tau_ = std::abs(tau_);
  }
  // 'tau' gets reset in case of non-convergence
  this->check_convergence(xv, yv);
}

template<>
inline void mcdata<VectorData>::find_conv_and_tau(void) const 
{
  throw std::runtime_error("mcdata<T>::find_conv_and_tau: this version not for 'VectorData' type");
}

template<typename T>
inline void mcdata<T>::find_conv_and_tau(const unsigned& n) const 
{ 
  throw std::runtime_error("mcdata<T>::find_conv_and_tau: this version is only for 'VectorData' type");
}

template<>
inline void mcdata<VectorData>::find_conv_and_tau(const unsigned& n) const 
{ 
  this->finalize();
  std::vector<double> xv, yv;
  unsigned level_n = dcorr_level_;
  if (dcorr_level_ == 2) level_n = dcorr_level_-2;
  else if (dcorr_level_ >= 3) level_n = dcorr_level_-3;
  for (unsigned i=level_n; i<=dcorr_level_; ++i) {
    xv.push_back(static_cast<double>(i));
    VectorData stddev = this->operator[](i).stddev();
    yv.push_back(static_cast<double>(stddev(n)));
  }
  // assuming convergence
  double stddev_0 = static_cast<double>(this->operator[](0).stddev()(n));
  double stddev_d = static_cast<double>(this->operator[](dcorr_level_).stddev()(n));
  if (stddev_0 < 1.0E-12) tau_ = -1.0;
  else {
    double r = stddev_d/stddev_0;
    tau_ = 0.5*( r*r - 1.0);
    tau_ = std::abs(tau_);
  }
  // 'tau' gets reset in case of non-convergence
  this->check_convergence(xv, yv);
  if (error_converged_ == "NOT_CONVD") tau_ = -1.0;
}

template<typename T>
inline void mcdata<T>::check_convergence(const std::vector<double>& xv, const std::vector<double>& yv) const 
{
  error_converged_ = "NOT_CONVD";
  if (xv.size() >= 3) {
    // slope of a least square fit straight line through these points
    double a1 = static_cast<double>(xv.size());
    double b1 = 0.0;
    for (const auto& x : xv) b1 += x;
    double c1 = 0.0;
    for (const auto& y : yv) c1 += y;
    double a2 = b1;
    double b2 = 0.0;
    for (const auto& x : xv) b2 += x*x;
    double c2 = 0.0;
    for (unsigned i=0; i<xv.size(); ++i) c2 += xv[i]*yv[i];
    double slope = (c2 - c1*a2/a1)/(b2 - b1*a2/a1);
    //if (std::abs(slope) < 1.0E-4) {
    if (std::abs(slope) < 0.1 * yv.back()) error_converged_ = "CONVERGED";
  }
} 

/*template<>
inline int mcdata<VectorData>::find_conv_and_tau(void) const 
{ 
  return 0;
} */
template<typename T>
inline std::string mcdata<T>::result_str(void) const 
{ 
  std::ostringstream os;
  os << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  os << std::setw(15) << this->mean() << std::fixed << std::setw(10) << this->stddev();
  return os.str();
}

template<>
inline std::string mcdata<VectorData>::result_str(void) const 
{ 
  throw std::runtime_error("mcdata<T>::result_str: this version is not for 'VectorData' types");
}

template<typename T>
inline std::string mcdata<T>::result_str(const int& n) const 
{ 
  throw std::runtime_error("mcdata<T>::result_str: this version is only for 'VectorData' types");
}

template<>
inline std::string mcdata<VectorData>::result_str(const int& n) const 
{ 
  std::ostringstream os;
  double mean, stddev;
  if (n<0) {
    mean = this->mean().sum();
    stddev = this->stddev().sum();
  }
  else {
    mean = this->mean()(n);
    stddev = this->stddev()(n);
  }
  os << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  os << std::setw(15) << mean << std::fixed << std::setw(10) << stddev;
  return os.str();
}

template<typename T>
inline std::string mcdata<T>::conv_str(const int& n) const 
{ 
  if (n < 0) this->find_conv_and_tau();
  else this->find_conv_and_tau(n);
  std::ostringstream os;
  os << std::scientific << std::uppercase << std::setprecision(6) << std::right;
  os << std::setw(10) << this->num_samples();
  os << std::setw(11) << error_converged_;
  os << std::setw(7) << std::setprecision(2) << std::right << std::fixed << tau_;
  os << std::resetiosflags(std::ios_base::floatfield) << std::nouppercase; 
  return os.str();
}


template<typename T>
inline void mcdata<T>::show_statistic(std::ostream& os) const
{
  auto bin = top_bin;
  os << name_ << " Statistic:\n";
  unsigned i = 0;
  while (bin->num_samples()>1) {
    os <<"bin-"<<std::setw(2)<<i++<<": "<< std::right<<std::setw(8)<<bin->num_samples(); 
    os << std::scientific << std::uppercase << std::setw(14) << bin->mean(); 
    os << std::setw(6) << "(+/-) " << std::fixed << std::setw(10) << bin->stddev(); 
    os << "\n";
    bin++;
  }
  os << std::right << std::resetiosflags(std::ios_base::floatfield) << std::nouppercase; 
}

template<>
inline void mcdata<VectorData>::show_statistic(std::ostream& os) const
{
  auto bin = top_bin;
  os << name_ << " Statistic:\n";
  unsigned i = 0;
  while (bin->num_samples()>1) {
    os <<"bin-"<<std::setw(2)<<i++<<": "<< std::right<<std::setw(8)<<bin->num_samples(); 
    VectorData mean = bin->mean();
    VectorData stddev = bin->stddev();
    for (int n=0; n<mean.rows(); ++n) {
      os << std::scientific << std::uppercase << std::setw(14) << mean(n); 
      os << std::setw(6) << "(+/-) " << std::fixed << std::setw(10) << stddev(n); 
    }
    os << "\n";
    bin++;
  }
  os << std::right << std::resetiosflags(std::ios_base::floatfield) << std::nouppercase; 
}


template<typename T>
inline std::ostream& operator<<(std::ostream& os, const mcdata<T>& obs)
{
  using namespace std;
  streamsize dp = cout.precision(); 
  os << scientific << uppercase << setprecision(6);
  os << obs.name() << " = " << setw(12) << obs.mean() << setw(8) << "(+/-) ";
  os << fixed << setw(8) << obs.stddev() << " (samples = " << obs.num_samples() << ")\n";
  if (obs.show_statistic_) {
    os << "\n";
    obs.show_statistic(os);
    obs.show_statistic_ = false;
  }
  os << resetiosflags(ios_base::floatfield) << nouppercase << setprecision(dp);
  return os;
}

template<>
inline std::ostream& operator<<(std::ostream& os, const mcdata<VectorData>& obs)
{
  using namespace std;
  streamsize dp = cout.precision(); 
  os << scientific << uppercase << setprecision(6);
  os << "#" << std::string(36, '-') << "\n";
  os << "# " << obs.name() << ":" << "(samples = " << obs.num_samples() << ")\n";
  os << left << setw(7) << "#i" << setw(16) << "mean" << setw(12) << "err" << "\n";
  os << "#" << std::string(36, '-') << "\n";
  auto mean = obs.mean();
  auto stddev = obs.stddev();
  for (int i=0; i<mean.rows(); ++i) {
    os << " ";
    os << scientific << setw(6) << i << setw(16) << mean(i); // << setw(6) << "(+/-)";
    os << fixed << setw(12) << stddev(i) << "\n";
  }
  if (obs.show_statistic_) {
    os << "\n";
    obs.show_statistic(os);
    obs.show_statistic_ = false;
  }
  os << right << resetiosflags(ios_base::floatfield) << nouppercase << setprecision(dp);
  return os;
}

using RealObservableData = mcdata<double>;
using VectorObservableData = mcdata<VectorData>;

} // end namespace mc

#endif

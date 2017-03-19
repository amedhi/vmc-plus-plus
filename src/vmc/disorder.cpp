/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-13 11:23:45
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-17 09:44:37
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
//#define BOOST_FILESYSTEM_NO_DEPRECATED
//#include <boost/filesystem.hpp>
#include <cstdlib>
#include "./disorder.h"

namespace vmc {

SiteDisorder::SiteDisorder(const input::Parameters& inputs)
{
  int nowarn;
  exists_ = inputs.set_value("put_site_disorder",false,nowarn);
}

int SiteDisorder::init(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model, RandomGenerator& rng)
{
  if (!model.have_disorder_term()) return 0;

  num_sites_ = graph.num_sites();
  num_configs_ = inputs.set_value("disorder_configs",1);
  bandwidth_ = inputs.set_value("disorder_bandwidth",1.0);
  disorder_pot_.resize(num_configs_);
  for (auto& v : disorder_pot_) v.resize(num_sites_);
  // disorder potential file path
  datafile_prefix_ = inputs.set_value("disorderinput_prefix","./disorder");
  if (datafile_prefix_.back()!='/') datafile_prefix_ += "/";
  datafile_prefix_ += model.signature_str();
  // load disorder potential from file, if available, otherwise 
  // generate the potentials
  std::string path = datafile_prefix_+"/cnf";
  for (int i=0; i<num_configs_; ++i) {
    std::ostringstream cnf_id;
    cnf_id.width(3);
    cnf_id.fill('0');
    cnf_id << i;
    std::string config_path(path+cnf_id.str()+"/");
    std::string fname(config_path+"site_disorder.txt");
    std::ifstream fin(fname);
    if (fin.is_open()) {
      double x;
      int n = 0;
      while (fin>>x && n<num_sites_) {
        disorder_pot_[i][n++] = x;
      }
      fin.close();
      if (n != num_sites_) {
        throw std::range_error("SiteDisorder::init: disorder data size mismatch");
      }
    }
    else {
      //boost::filesystem::path dir(config_path);
      //boost::filesystem::create_directory(dir);
      std::string cmd("mkdir -p "+config_path);
      std::system(cmd.c_str());
      std::ofstream fout(fname);
      //std::cout << fname << "\n";
      if (!fout.is_open()) 
        throw std::runtime_error("SiteDisorder::init: file open failed");
      // random site potential
      fout<<std::scientific<<std::setprecision(12)<<std::uppercase<<std::right;
      for (unsigned n=0; n<num_sites_; ++n) {
        double V = bandwidth_ * (rng.random_real()-0.50);
        //fout << std::setw(6) << n << std::setw(22) << V << std::endl;
        fout << std::setw(22) << V << std::endl;
      }
      fout.close();
    }
    //config += path.str()+"/";
  }
  std::cout << "SiteDisorder\n"; getchar();

  return 0;
}
  




} // end namespace vmc


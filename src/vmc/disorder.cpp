/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-13 11:23:45
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-14 00:30:08
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
//#define BOOST_FILESYSTEM_NO_DEPRECATED
//#include <boost/filesystem.hpp>
#include <cstdlib>
#include "./disorder.h"

namespace vmc {

int SiteDisorder::init(const input::Parameters& inputs, const lattice::LatticeGraph& graph,
    const model::Hamiltonian& model)
{
  if (!model.have_disorder_term()) return 0;

  num_sites_ = graph.num_sites();
  num_configs_ = inputs.set_value("disorder_configs",1);
  bandwidth_ = inputs.set_value("disorder_bandwidth",1.0);
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
    std::string fname(config_path+"disorder_pot.txt");
    std::ifstream fin(fname);
    if (fin.is_open()) {
    }
    else {
      //boost::filesystem::path dir(config_path);
      //boost::filesystem::create_directory(dir);
      std::string cmd("mkdir -p "+config_path);
      std::system(cmd.c_str());
      std::ofstream fout(fname);
      std::cout << fname << "\n";
      if (!fout.is_open()) 
        throw std::runtime_error("SiteDisorder::input: file open failed");
      for (unsigned n=0; n<num_sites_; ++n)
        fout << bandwidth_ << std::endl;
      fout.close();
    }
    //config += path.str()+"/";
  }
  std::cout << "SiteDisorder\n"; getchar();

  return 0;
}
  




} // end namespace vmc


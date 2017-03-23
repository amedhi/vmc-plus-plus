/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-22 22:46:55
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-24 00:15:14
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./disordered_sc.h"

namespace var {

DisorderedSC::DisorderedSC(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
  : GroundState(true)
{
  init(inputs, graph);
}

int DisorderedSC::init(const input::Parameters& inputs, const lattice::LatticeGraph& graph)
{
  // sites & bonds
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  // particle number
  set_particle_num(inputs);
  // infinity limit
  large_number_ = 1.0E+2;

  // the hamiltonian to be stored in Sparse matrix formats
  quadratic_ham_.resize(num_sites_,num_sites_);
  pairing_ham_.resize(num_sites_,num_sites_);
  quadratic_coeffs_.resize(num_sites_+num_bonds_);
  pairing_coeffs_.resize(num_bonds_);
  delta_.resize(num_sites_,num_sites_);
  work_.resize(num_sites_,num_sites_);
  phi_k_.resize(num_sites_);

  // hopping parameters
  // variational parameters
  varparms_.clear();
  std::string name;
  double defval, lb, ub;
  // site potentials
  mu_start_ = varparms_.size();
  for (unsigned i=0; i<num_sites_; ++i) {
    name = "mu_" + std::to_string(i);
    varparms_.add(name, defval=0.0, lb=-4.0, ub=4.0);
  }
  t_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "t_" + std::to_string(i);
    varparms_.add(name, defval=-1.0, lb=-2.0, ub=2.0);
  }
  delta_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "delta_" + std::to_string(i);
    varparms_.add(name, defval=0.5, lb=-4.0, ub=4.0);
  }
  num_varparms_ = varparms_.size();

  // Hamiltonian coefficients
  // mu-values
  for (unsigned i=0; i<num_sites_; ++i) {
    quadratic_coeffs_[i] = EigenTriplet(i,i,varparms_[mu_start_+i].value());
  }
  // t & delta values (only lower triangular elements)
  unsigned i = num_sites_;
  unsigned j = 0;
  unsigned it = t_start_;
  unsigned id = delta_start_;
  for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int site_i = graph.source(b);
    int site_j = graph.target(b);
    // t-values
    double t = varparms_[it++].value();
    double delta = varparms_[id++].value();
    if (site_i >= site_j) {
      quadratic_coeffs_[i++] = EigenTriplet(site_i,site_j,t);
      pairing_coeffs_[j++] = EigenTriplet(site_i,site_j,delta);
    }
    else {
      quadratic_coeffs_[i++] = EigenTriplet(site_j,site_i,t);
      pairing_coeffs_[j++] = EigenTriplet(site_j,site_i,delta);
    }
  }
  /*for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int site_i = graph.source(b);
    int site_j = graph.target(b);
    double t = varparms_[it++].value();
    quadratic_coeffs_[i++] = EigenTriplet(site_i,site_j,t);
    quadratic_coeffs_[i++] = EigenTriplet(site_j,site_i,t);
    double delta = varparms_[id++].value();
    pairing_coeffs_[j++] = EigenTriplet(site_i,site_j,delta);
    pairing_coeffs_[j++] = EigenTriplet(site_j,site_i,delta);
  }*/
  // set the sparse matrices
  quadratic_ham_.setFromTriplets(quadratic_coeffs_.begin(), quadratic_coeffs_.end());
  pairing_ham_.setFromTriplets(pairing_coeffs_.begin(), pairing_coeffs_.end());
  // check whether sparse matrix diagonalization is okay
  /*delta_.setZero();
  for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int i = graph.source(b);
    int j = graph.target(b);
    delta_(i,j) = -1.0;
    delta_(j,i) = -1.0;
  }
  es_quad_.compute(quadratic_ham_);
  es_pair_.compute(delta_);
  std::cout << es_quad_.eigenvectors() - es_pair_.eigenvectors() << "\n\n";
  getchar();*/


  return 0;
}

void DisorderedSC::update(const input::Parameters& inputs)
{
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
}

void DisorderedSC::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  // update mu values
  unsigned i = 0;
  for (auto& coeff : quadratic_coeffs_) {
    coeff = EigenTriplet(coeff.row(), coeff.col(), pvector[start_pos+i]);
    ++i;
  }
  for (auto& coeff : pairing_coeffs_) {
    coeff = EigenTriplet(coeff.row(), coeff.col(), pvector[start_pos+i]);
    ++i;
  }
  // new Hamiltonian matrices
  quadratic_ham_.setFromTriplets(quadratic_coeffs_.begin(), quadratic_coeffs_.end());
  pairing_ham_.setFromTriplets(pairing_coeffs_.begin(), pairing_coeffs_.end());
}

void DisorderedSC::get_wf_amplitudes(Matrix& psi) 
{
  // diagonalize the quadratic part
  es_quad_.compute(quadratic_ham_);
  // work = \Delta + Transpose(Delta)
  work_ = pairing_ham_ + Eigen::SparseMatrix<double>(pairing_ham_.transpose());
  // transform the pairing part
  delta_ = es_quad_.eigenvectors().adjoint() * work_ * 
      es_quad_.eigenvectors().conjugate();
  //std::cout << work_ << "\n"; getchar();
  // for INTRABAND pairing, phi_a = v_a/u_a for each band a
  for (unsigned k=0; k<num_sites_; ++k) {
    double ek = es_quad_.eigenvalues()[k];
    double delta_k = delta_(k,k);
    double deltak_sq = delta_k * delta_k;
    if (deltak_sq<1.0E-12 && ek>0.0) {
      phi_k_[k] = large_number_;
    }
    else {
      phi_k_[k] = delta_k/(-ek+std::sqrt(ek*ek+deltak_sq));
    }
    //std::cout << ek << "\n";
    //std::cout << delta_k << "\n";
    //std::cout << phi_k_[k] << "\n";
    //getchar();
  }
  // bcs ampitudes in site basis (delta used as work array)
  for (unsigned i=0; i<num_sites_; ++i) 
    delta_.col(i) = es_quad_.eigenvectors().col(i) * phi_k_[i];
  psi = delta_ * es_quad_.eigenvectors().transpose();
}

void DisorderedSC::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  for (unsigned i=0; i<t_start_; ++i) {
    double h = varparms_[i].diff_h();
    double inv_2h = 0.5/h;
    quadratic_ham_.coeffRef(i,i) += h;
  }
}

void DisorderedSC::get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
}

void DisorderedSC::get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k)
{
}



} // end namespace var
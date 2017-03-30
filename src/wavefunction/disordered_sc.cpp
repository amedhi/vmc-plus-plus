/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-22 22:46:55
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-30 23:07:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./disordered_sc.h"
#include <fstream>

//#define HACK_ON

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
  large_number_ = 1.0E+6;

  // the hamiltonian to be stored in Sparse matrix formats
  quadratic_ham_.resize(num_sites_,num_sites_);
  pairing_ham_.resize(num_sites_,num_sites_);
  quadratic_coeffs_.resize(num_sites_+num_bonds_);
  pairing_coeffs_.resize(num_bonds_);
  delta_.resize(num_sites_,num_sites_);
  phi_k_.resize(num_sites_);
  work_.resize(num_sites_,num_sites_);
  psi_work_.resize(num_sites_,num_sites_);
  //psi_work2_.resize(num_sites_,num_sites_);

  // variational parameters
  varparms_.clear();
  std::string name;
  double defval, lb, ub, h;
  // site potentials
  mu_start_ = varparms_.size();
  for (unsigned i=0; i<num_sites_; ++i) {
    name = "mu_" + std::to_string(i);
    varparms_.add(name, defval=-0.50, lb=-1.5, ub=1.5, h=0.2);
  }
  t_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "t_" + std::to_string(i);
    varparms_.add(name, defval=1.0, lb=0.5, ub=1.2, h=0.01);
  }
  delta_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "delta_" + std::to_string(i);
    varparms_.add(name, defval=0.4, lb=0.0, ub=0.8, h=0.01);
  }
  num_varparms_ = varparms_.size();

  // Hamiltonian coefficients
  // mu-values
  for (unsigned i=0; i<num_sites_; ++i) {
    double mu = varparms_[mu_start_+i].value();
    quadratic_coeffs_[i] = MatrixElem(i,i,mu);
  }
  // t & delta values (only lower triangular elements)
  unsigned i = num_sites_;
  unsigned j = 0;
  unsigned it = t_start_;
  unsigned id = delta_start_;
  //std::cout << "** DisorderedSC::init: ---TEMPORARY HACK-----\n";
  for (auto b=graph.bonds_begin(); b!=graph.bonds_end(); ++b) {
    int site_i = graph.source(b);
    int site_j = graph.target(b);
    // t-values
    double t = varparms_[it++].value();
    if (site_i >= site_j) 
      quadratic_coeffs_[i++] = MatrixElem(site_i,site_j,t,graph.bond_sign(b));
    else 
      quadratic_coeffs_[i++] = MatrixElem(site_j,site_i,t,graph.bond_sign(b));
    // delta values (assuming DWAVE SC order)
    double delta = varparms_[id++].value();
    switch (graph.bond_type(b)) {
      case 0:
        pairing_coeffs_[j++] = MatrixElem(site_i,site_j,delta,+1);
        break;
      case 1:
        pairing_coeffs_[j++] = MatrixElem(site_i,site_j,delta,-1);
        break;
      default:
        pairing_coeffs_[j++] = MatrixElem(site_i,site_j,delta,+1);
        break;
    }
  }
  // quadratic hamilonian (only LOWER TRIANGULAR part)
  quadratic_ham_.setZero();
  for (auto& coeff : quadratic_coeffs_) {
    quadratic_ham_(coeff.row(), coeff.col()) = -coeff.value()*coeff.bond_phase();
  }
  // pairing Hamiltonian
  pairing_ham_.setZero();
  for (auto& coeff : pairing_coeffs_) {
    double val = -coeff.value()*coeff.bond_phase();
    pairing_ham_(coeff.row(), coeff.col()) = val;
    pairing_ham_(coeff.col(), coeff.row()) = val;
  }

  // HACK
#ifdef HACK_ON
  num_varparms_ = 1;
  varparms_.resize(1);
#endif

  return 0;
}

void DisorderedSC::update(const input::Parameters& inputs)
{
  // HACK
#ifdef HACK_ON
  double x = inputs.set_value("varp", 1.0);
  for (unsigned i=0; i<num_sites_; i++) {
    //varparms_[i].change_value(x);
    auto coeff = quadratic_coeffs_[i];
    coeff.change_value(x);
    quadratic_ham_(coeff.row(),coeff.col()) = -coeff.value();
    //std::cout << coeff.row() << " " << coeff.col() << " " << -x << "\n";
  }
  //std::cout << quadratic_ham_.diagonal().transpose() << "\n"; getchar();
#endif
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
}

void DisorderedSC::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  // update matrix elements
  // '-' sign because H = -\sum t_{ij} cicj - mu\sum_ni -sum_{ij}Delta_{ij}cicj etc
  unsigned i = 0;
  for (auto& coeff : quadratic_coeffs_) {
    //double x = pvector[start_pos+i];
    //varparms_[i].change_value(x);
    coeff.change_value(pvector[start_pos+i]);
    double val = -coeff.value()*coeff.bond_phase();
    quadratic_ham_(coeff.row(), coeff.col()) = val;
    ++i;
  }
  for (auto& coeff : pairing_coeffs_) {
    //double delta = pvector[start_pos+i];
    //varparms_[i].change_value(delta);
    coeff.change_value(pvector[start_pos+i]);
    double val = -coeff.value()*coeff.bond_phase();
    pairing_ham_(coeff.row(), coeff.col()) = val;
    pairing_ham_(coeff.col(), coeff.row()) = val;
    ++i;
  }
}

void DisorderedSC::get_wf_amplitudes(Matrix& psi) 
{
  // diagonalize the quadratic part
  es_quad_.compute(quadratic_ham_);
  // transform the pairing part
  delta_ = es_quad_.eigenvectors().adjoint() * pairing_ham_ * 
      es_quad_.eigenvectors().conjugate();
  // wf amplitudes
  get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
    delta_.diagonal(), psi);
  //std::cout << psi  << "\n"; getchar();
  //psi_work_ = psi;
  //std::cout << es_quad_.eigenvalues() << "\n"; getchar();
  //std::cout << delta_.diagonal() << "\n"; getchar();
}

void DisorderedSC::hack_gradient(std::vector<Matrix>& psi_gradient) 
{
  //ComplexMatrix tmp(num_sites_,num_sites_);
  unsigned i = 0;
  double h = 0.2; // 5.0E-1;// varparms_[0].diff_h();
  double inv_2h = 0.5/h;
  auto coeff = quadratic_coeffs_[i];
  double x = coeff.value();
  std::cout << "x = " << x << "\n";
  std::cout << "h = " << h << "\n";
  std::cout << "1/2h = " << inv_2h << "\n"; 
  for (unsigned j=0; j<num_sites_; ++j)
  quadratic_ham_(j,j) = -(x+h);
  //quadratic_ham_(coeff.col(),coeff.row()) = -(x+h);
  //std::cout << quadratic_ham_ << "\n"; getchar();
  get_wf_amplitudes(psi_gradient[i]);
  for (unsigned j=0; j<num_sites_; ++j)
  quadratic_ham_(j,j) = -(x-h);
  //quadratic_ham_(coeff.col(),coeff.row()) = -(x-h);
  //std::cout << quadratic_ham_ << "\n"; getchar();
  get_wf_amplitudes(psi_work_);
  // restore the hamtonian
  for (unsigned j=0; j<num_sites_; ++j)
  quadratic_ham_(j,j) = -x;
  //quadratic_ham_(coeff.col(),coeff.row()) = -x;
  // gradient
  psi_gradient[i] -= psi_work_;
  //std::cout << psi_gradient[i] << "\n"; getchar();
  //std::cout << "-----hi-------" << "\n"; getchar();
  psi_gradient[i] *= inv_2h;
}

void DisorderedSC::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{

#ifdef HACK_ON
  return hack_gradient(psi_gradient);
  std::cout << "---NO GO HERE-------" << "\n";
#endif

  // derivative wrt quadratic parameters
  unsigned i = 0;
  for (auto& coeff : quadratic_coeffs_) {
    double val = quadratic_ham_(coeff.row(),coeff.col());
    double x = coeff.value();
    double h = varparms_[i].diff_h();
    double inv_2h = 0.5/h;
    quadratic_ham_(coeff.row(),coeff.col()) = -(x+h)*coeff.bond_phase();
    get_wf_amplitudes(psi_gradient[i]);
    quadratic_ham_(coeff.row(),coeff.col()) = -(x-h)*coeff.bond_phase();
    get_wf_amplitudes(psi_work_);
    // restore the hamtonian
    quadratic_ham_(coeff.row(),coeff.col()) = val;
    // gradient
    psi_gradient[i] -= psi_work_;
    psi_gradient[i] *= inv_2h;
    ++i;
  }
  // derivative wrt pairing parameters
  // quadratic part remains constant
  es_quad_.compute(quadratic_ham_);
  // paring matrix elements
  for (auto& coeff : pairing_coeffs_) {
    double val = pairing_ham_(coeff.row(),coeff.col());
    double x = coeff.value();
    double h = varparms_[i].diff_h();
    double inv_2h = 0.5/h;
    // wavefunction 'at (x+h)'
    double val_plus = -(x+h) * coeff.bond_phase();
    pairing_ham_(coeff.row(),coeff.col()) = val_plus;
    pairing_ham_(coeff.col(),coeff.row()) = val_plus;
    delta_ = es_quad_.eigenvectors().adjoint() * pairing_ham_ * 
             es_quad_.eigenvectors().conjugate();
    get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
      delta_.diagonal(), psi_gradient[i]);
    // wavefunction 'at (x-h)'
    double val_minus = -(x-h) * coeff.bond_phase();
    pairing_ham_(coeff.row(),coeff.col()) = val_minus;
    pairing_ham_(coeff.col(),coeff.row()) = val_minus;
    delta_ = es_quad_.eigenvectors().adjoint() * pairing_ham_ * 
             es_quad_.eigenvectors().conjugate();
    get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
      delta_.diagonal(), psi_work_);
    // restore the elements
    pairing_ham_(coeff.row(),coeff.col()) = val;
    pairing_ham_(coeff.col(),coeff.row()) = val;
    // gradient
    psi_gradient[i] -= psi_work_;
    psi_gradient[i] *= inv_2h;
    ++i;
  }
}

void DisorderedSC::get_pair_amplitudes_sitebasis(const Eigen::VectorXd& es_eigenvalues, 
  const RealMatrix& es_eigenvectors, const Eigen::VectorXd& delta, Matrix& psi)
{
  // for INTRABAND pairing, phi_a = v_a/u_a for each band a
  //std::cout << delta.diagonal().transpose() << "\n";
  for (unsigned k=0; k<num_sites_; ++k) {
    double ek = es_eigenvalues[k];
    auto delta_k = delta(k);
    if (std::abs(delta_k)<1.0E-12 && ek<0.0) {
      if (delta_k>=0.0) phi_k_[k] = large_number_;
      else phi_k_[k] = -large_number_;
      //std::cout << "delta_k = " << delta_k << "\n";
      //std::cout << "ek = " << ek << "\n";
      //std::cout << "large_number_\n";
      //getchar();
    }
    /*
    else if (std::abs(delta_k)<1.0E-12 && std::abs(ek)<1.0E-12) {
      phi_k_[k] = 1.0;
    }*/
    else {
      phi_k_[k] = delta_k/(ek+std::sqrt(ek*ek+delta_k*delta_k));
    }
    //std::cout << ek << "\n";
    //std::cout << delta_k << "\n";
    //std::cout << phi_k_[k] << "\n";
    //getchar();
  }
  //getchar();
  // pair ampitudes in site basis 
  for (unsigned i=0; i<num_sites_; ++i) {
    work_.col(i) = es_eigenvectors.col(i) * phi_k_[i];
  }
  work_ *= es_eigenvectors.transpose();
  double one_by_nk = 1.0/static_cast<double>(num_sites_);
  for (unsigned i=0; i<num_sites_; ++i) {
    for (unsigned j=0; j<num_sites_; ++j) {
      psi(i,j) = ampl_part(work_(i,j)) * one_by_nk;
      //std::cout << psi_work2_(i,j) << "\n"; getchar();
    }
  }
}



} // end namespace var
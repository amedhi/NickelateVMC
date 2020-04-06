/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-22 22:46:55
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-22 23:35:21
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./disordered_sc.h"
#include <fstream>

#define HACK_ON
//#define CHECK_T
#define CHECK_D

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
  mu_start_ = 0;
  for (unsigned i=0; i<num_sites_; ++i) {
    name = "mu_" + std::to_string(i);
    //varparms_.add(name, defval=-0.50, lb=-2.5, ub=2.5, h=0.2);
    varparms_.add(name, defval=-0.50, lb=-2.5, ub=2.5, h=0.1);
  }
  t_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "t_" + std::to_string(i);
    varparms_.add(name, defval=1.0, lb=0.1, ub=5.0, h=0.1);
  }
  delta_start_ = varparms_.size();
  for (unsigned i=0; i<num_bonds_; ++i) {
    name = "delta_" + std::to_string(i);
    //varparms_.add(name, defval=0.4, lb=0.0, ub=1.5, h=0.1);
    varparms_.add(name, defval=1.0, lb=0.00, ub=2.0, h=0.1);
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
#ifdef CHECK_T
  varparms_.clear();
  num_varparms_ = 1;
  double x = inputs.set_value("varp", 1.0);
  varparms_.add("tv", defval=x, lb=0.1, ub=10.5, h=0.1);
  for (unsigned i=t_start_; i<delta_start_; i++) {
  //for (unsigned i=t_start_; i<t_start_+1; i++) {
    auto coeff = quadratic_coeffs_[i];
    coeff.change_value(varparms_[0].value());
    quadratic_ham_(coeff.row(),coeff.col()) = -coeff.value()*coeff.bond_phase();
  }
#endif
#ifdef CHECK_D
  varparms_.clear();
  num_varparms_ = 1;
  double x = inputs.set_value("varp", 1.0);
  varparms_.add("d", defval=x, lb=0.0, ub=2.0, h=0.1);
  //for (auto& coeff : pairing_coeffs_) coeff.change_value(x);
  pairing_coeffs_[0].change_value(x);
  //for (int i=0; i<10; ++i) pairing_coeffs_[i].change_value(x);
#endif

  return 0;
}

void DisorderedSC::update(const input::Parameters& inputs)
{
  // HACK
#ifdef CHECK_T
  unsigned start_p = t_start_;
  unsigned end_p = delta_start_;
  double x = inputs.set_value("varp", 1.0);
  for (unsigned i=start_p; i<end_p; i++) {
  //for (unsigned i=start_p; i<start_p+1; i++) {
    //varparms_[i].change_value(x);
    auto coeff = quadratic_coeffs_[i];
    coeff.change_value(x);
    quadratic_ham_(coeff.row(),coeff.col()) = -coeff.value()*coeff.bond_phase();
    //std::cout << coeff.row() << " " << coeff.col() << " " << -x << "\n";
  }
  //std::cout << quadratic_ham_.diagonal().transpose() << "\n"; getchar();
  //std::cout << quadratic_ham_ << "\n"; getchar();
#endif
#ifdef CHECK_D
  double x = inputs.set_value("varp", 1.0);
  //for (auto& coeff : pairing_coeffs_) coeff.change_value(x);
  pairing_coeffs_[0].change_value(x);
  //for (int i=0; i<10; ++i) pairing_coeffs_[i].change_value(x);
#endif

  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
}

void DisorderedSC::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  // HACK
#ifdef CHECK_T
  unsigned start_p = t_start_;
  unsigned end_p = delta_start_;
  double x = pvector[start_pos];
  //for (unsigned i=0; i<num_sites_; i++) {
  for (unsigned i=start_p; i<start_p+1; i++) {
    //varparms_[i].change_value(x);
    auto coeff = quadratic_coeffs_[i];
    coeff.change_value(x);
    quadratic_ham_(coeff.row(),coeff.col()) = -coeff.value()*coeff.bond_phase();
    //std::cout << coeff.row() << " " << coeff.col() << " " << -x << "\n";
  }
  //std::cout << quadratic_ham_.diagonal().transpose() << "\n"; getchar();
  //std::cout << quadratic_ham_ << "\n"; getchar();
  return;
#endif
#ifdef CHECK_D
  double x = pvector[start_pos];
  pairing_coeffs_[0].change_value(x);
  return;
#endif


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
  //std::cout << "\n" << pairing_ham_ << "\n"; getchar();
}

void DisorderedSC::get_wf_amplitudes(Matrix& psi) 
{
  // diagonalize the quadratic part
  es_quad_.compute(quadratic_ham_);

  // transform the pairing part
  //delta_ = es_quad_.eigenvectors().adjoint() * pairing_ham_ * 
  //    es_quad_.eigenvectors().conjugate();
  //std::cout << "\n" << quadratic_ham_ << "\n"; getchar();

  // delta_ calculation above was wrong?
  for (unsigned k=0; k<num_sites_; ++k) {
    double sum = 0.0;
    for (auto& coeff : pairing_coeffs_) {
      int i = coeff.row();
      int j = coeff.col();
      double delta_ij = coeff.value()*coeff.bond_phase();
      sum += delta_ij * es_quad_.eigenvectors()(i,k) * es_quad_.eigenvectors()(j,k); 
    }
    delta_(k,k) = -sum;
  }

  // wf amplitudes
  get_pair_amplitudes_sitebasis(es_quad_.eigenvalues(), es_quad_.eigenvectors(),
    delta_.diagonal(), psi);
  //std::cout << psi  << "\n"; getchar();
  //psi_work_ = psi;
  //std::cout << es_quad_.eigenvalues() << "\n"; getchar();
  //std::cout << delta_.diagonal().sum() << "\n"; getchar();
}

void DisorderedSC::hack_gradient(std::vector<Matrix>& psi_gradient) 
{
#ifdef CHECK_T
  //ComplexMatrix tmp(num_sites_,num_sites_);
  unsigned start_p = t_start_;
  unsigned end_p = delta_start_;
  double h = varparms_[0].diff_h();
  double inv_2h = 0.5/h;
  auto coeff = quadratic_coeffs_[start_p];
  double save = quadratic_ham_(coeff.row(),coeff.col());
  //for (unsigned i=start_p; i<end_p; i++) {
  for (unsigned i=start_p; i<start_p+1; i++) {
    auto coeff = quadratic_coeffs_[i];
    double x = coeff.value();
    quadratic_ham_(coeff.row(),coeff.col()) = -(x+h)*coeff.bond_phase();
  }
  get_wf_amplitudes(psi_gradient[0]);
  //for (unsigned i=start_p; i<end_p; i++) {
  for (unsigned i=start_p; i<start_p+1; i++) {
    auto coeff = quadratic_coeffs_[i];
    double x = coeff.value();
    quadratic_ham_(coeff.row(),coeff.col()) = -(x-h)*coeff.bond_phase();
  }
  get_wf_amplitudes(psi_work_);
  // restore the hamtonian
  //for (unsigned i=start_p; i<end_p; i++) {
  for (unsigned i=start_p; i<start_p+1; i++) {
    auto coeff = quadratic_coeffs_[i];
    quadratic_ham_(coeff.row(),coeff.col()) = save;
  }
  // gradient
  psi_gradient[0] -= psi_work_;
  psi_gradient[0] *= inv_2h;
  //std::cout << "-----hi-------" << "\n"; getchar();
#endif


#ifdef CHECK_D
  // numerical gradient
/*
  double h = varparms_[0].diff_h();
  double inv_2h = 0.5/h;
  double x = pairing_coeffs_[0].value();
  //for (auto& coeff : pairing_coeffs_) coeff.change_value(x+h);
  pairing_coeffs_[0].change_value(x+h);
  get_wf_amplitudes(psi_gradient[0]);
  //for (auto& coeff : pairing_coeffs_) coeff.change_value(x-h);
  pairing_coeffs_[0].change_value(x-h);
  get_wf_amplitudes(psi_work_);
  //for (auto& coeff : pairing_coeffs_) coeff.change_value(x);
  pairing_coeffs_[0].change_value(x);
  // gradient
  psi_gradient[0] -= psi_work_;
  psi_gradient[0] *= inv_2h;
*/
  // analytical gradient
  es_quad_.compute(quadratic_ham_);
  // delta_k
  for (unsigned k=0; k<num_sites_; ++k) {
    double sum = 0.0;
    for (auto& coeff : pairing_coeffs_) {
      int i = coeff.row();
      int j = coeff.col();
      double delta_ij = coeff.value()*coeff.bond_phase();
      sum += delta_ij * es_quad_.eigenvectors()(i,k) * es_quad_.eigenvectors()(j,k); 
    }
    delta_(k,k) = -sum;
  }
  std::vector<double> alpha(num_sites_);
  for (unsigned k=0; k<num_sites_; ++k) {
    double ek = es_quad_.eigenvalues()[k];
    double delta_k = delta_(k,k);
    double Ek = std::sqrt(ek*ek+delta_k*delta_k);
    if (std::abs(delta_k)<1.0e-15 && ek<0.0) {
      alpha[k] = -large_number_;
      std::cout << "large_number_ = " << -large_number_; getchar();
    }
    else {
      alpha[k] = ek/(Ek*(ek + Ek));
    }
  }
  for (unsigned p=0; p<1; ++p) {
    int m = pairing_coeffs_[p].row();
    int n = pairing_coeffs_[p].col();
    for (unsigned i=0; i<num_sites_; ++i) {
      for (unsigned j=0; j<num_sites_; ++j) {
        double sum = 0.0;
        for (unsigned k=0; k<num_sites_; ++k) {
          sum += alpha[k] * es_quad_.eigenvectors()(m,k) * es_quad_.eigenvectors()(n,k) 
            * es_quad_.eigenvectors()(i,k) * es_quad_.eigenvectors()(j,k); 
        }
        psi_gradient[p](i,j) = sum;
        //std::cout << sum << "\n"; getchar();
      }
    }
  }
#endif
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
  // analytical gradient
  es_quad_.compute(quadratic_ham_);
  // delta_k
  for (unsigned k=0; k<num_sites_; ++k) {
    double sum = 0.0;
    for (auto& coeff : pairing_coeffs_) {
      int i = coeff.row();
      int j = coeff.col();
      double delta_ij = coeff.value()*coeff.bond_phase();
      sum += delta_ij * es_quad_.eigenvectors()(i,k) * es_quad_.eigenvectors()(j,k); 
    }
    delta_(k,k) = -sum;
  }
  std::vector<double> alpha(num_sites_);
  for (unsigned k=0; k<num_sites_; ++k) {
    double ek = es_quad_.eigenvalues()[k];
    double delta_k = delta_(k,k);
    double Ek = std::sqrt(ek*ek+delta_k*delta_k);
    if (std::abs(delta_k)<1.0e-15 && ek<0.0) {
      alpha[k] = -large_number_;
      //getchar();
    }
    else {
      alpha[k] = ek/(Ek*(ek + Ek));
    }
    //std::cout << "phi_k = " << phi_k_[k] << "\n";
    //getchar();
  }
  unsigned p = delta_start_;
  for (auto& coeff : pairing_coeffs_) {
    int m = coeff.row();
    int n = coeff.col();
    for (unsigned i=0; i<num_sites_; ++i) {
      for (unsigned j=0; j<num_sites_; ++j) {
        double sum = 0.0;
        for (unsigned k=0; k<num_sites_; ++k) {
          sum += alpha[k] * es_quad_.eigenvectors()(m,k) * es_quad_.eigenvectors()(n,k) 
            * es_quad_.eigenvectors()(i,k) * es_quad_.eigenvectors()(j,k); 
        }
        psi_gradient[p](i,j) = sum;
      }
    }
#ifdef CHECK_D
    return;
#endif
    ++p;
  }

  /*
  // numerical gradient
  for (auto& coeff : pairing_coeffs_) {
    double x = coeff.value();
    double h = varparms_[i].diff_h();
    double inv_2h = 0.5/h;
    coeff.change_value(x+h);
    get_wf_amplitudes(psi_gradient[i]);
    coeff.change_value(x-h);
    get_wf_amplitudes(psi_work_);
    // restore the elements
    coeff.change_value(x);
    // gradient
    psi_gradient[i] -= psi_work_;
    psi_gradient[i] *= inv_2h;
    //std::cout << psi_gradient[i] << "\n"; getchar();
    //std::cout << i << "\n"; getchar();
    ++i;
  }*/
  return;
  std::cout << "---NO GO HERE-------" << "\n";

  /*
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
  }*/
}

void DisorderedSC::get_gradient_pairing_coeff(std::vector<Matrix>& psi_gradient,
  const unsigned& start_pos)
{

} 

void DisorderedSC::get_pair_amplitudes_sitebasis(const Eigen::VectorXd& es_eigenvalues, 
  const RealMatrix& es_eigenvectors, const Eigen::VectorXd& delta, Matrix& psi)
{
  // for INTRABAND pairing, phi_a = v_a/u_a for each band a
  //std::cout << delta.diagonal().transpose() << "\n";
  for (unsigned k=0; k<num_sites_; ++k) {
    double ek = es_eigenvalues[k];
    auto delta_k = delta(k);
    if (std::abs(delta_k)<1.0e-12 && ek<0.0) {
      if (delta_k>=0.0) phi_k_[k] = large_number_;
      else phi_k_[k] = -large_number_;
      //std::cout << "delta_k = " << delta_k << "\n";
      //std::cout << "ek = " << ek << "\n";
      //std::cout << "large_number_\n";
      //getchar();
    }
    /*
    else if (std::abs(delta_k)<1.0e-12 && std::abs(ek)<1.0e-12) {
      phi_k_[k] = 1.0;
    }*/
    else {
      phi_k_[k] = delta_k/(ek+std::sqrt(ek*ek+delta_k*delta_k));
    }
    //std::cout << "k = " << k << "\n";
    //std::cout << "ek = " << ek << "\n";
    //std::cout << "delta_k = " << delta_k << "\n";
    //std::cout << "phi_k = " << phi_k_[k] << "\n";
    //getchar();
  }
  //getchar();
  // pair ampitudes in site basis 
  for (unsigned i=0; i<num_sites_; ++i) {
    for (unsigned j=0; j<num_sites_; ++j) {
      double sum = 0.0;
      for (unsigned k=0; k<num_sites_; ++k) {
        sum += phi_k_[k] * es_eigenvectors(i,k) * es_eigenvectors(j,k); 
      }
      psi(i,j) = sum;
    }
  }

  /*
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
  */
}



} // end namespace var
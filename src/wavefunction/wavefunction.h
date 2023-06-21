/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-11 10:28:23
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef WAVEFUNCTION_H
#define WAVEFUNCTION_H

#include <iostream>
#include <sstream>
#include <string>
#include <complex>
#include <vector>
//#include <map>
#include <memory>
#include <stdexcept>
#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "./mf_model.h"
#include "./matrix.h"
#include "./fermisea.h"
#include "./bcs_state.h"
#include "./disordered_sc.h"

namespace var {

enum class wf_type {normal, bcs_oneband, bcs_multiband, bcs_disordered};

// Wavefunction in 'Wannier space' (amplitudes in Wannier representation)
class Wavefunction 
{
public:
  //Wavefunction() {}
  Wavefunction(const lattice::Lattice& lattice, const input::Parameters& inputs,
    const model::Hamiltonian& model, const bool& site_disorder=false);
  ~Wavefunction() {}
  const MF_Order::pairing_t& pair_symmetry(void) const { return groundstate_->pair_symm(); }
  const MF_Order& mf_order(void) const { return *groundstate_; }
  const VariationalParms& varparms(void) const { return groundstate_->varparms(); }
  std::string info_str(void) const { return groundstate_->info_str(); } 
  int compute(const lattice::Lattice& lattice, const input::Parameters& inputs, 
    const bool& psi_gradient=false);
  int compute(const lattice::Lattice& lattice, const var::parm_vector& pvector,
    const unsigned& start_pos, const bool& psi_gradient=false);
  int recompute(const lattice::Lattice& lattice);
  //int compute_gradients(const lattice::Lattice& lattice);
  const int& num_upspins(void) const { return groundstate_->num_upspins(); }
  const int& num_dnspins(void) const { return groundstate_->num_dnspins(); }
  const double& hole_doping(void) const { return groundstate_->hole_doping(); }
  const basis::BlochBasis blochbasis(void) const { return groundstate_->blochbasis(); }
  std::string signature_str(void) const; 
  void get_vparm_names(std::vector<std::string>& names, unsigned start_pos) const; 
  void get_vparm_values(var::parm_vector& values, unsigned start_pos);
  void get_vparm_vector(std::vector<double>& vparm_values, unsigned start_pos);
  void get_vparm_lbound(var::parm_vector& lbounds, unsigned start_pos) const; 
  void get_vparm_ubound(var::parm_vector& ubounds, unsigned start_pos) const; 
  void get_amplitudes(Matrix& psi, const std::vector<int>& row,  
    const std::vector<int>& col) const;
  void get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const;
  void get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const;
  void get_amplitudes(amplitude_t& elem, const int& irow, const int& jcol) const;
  void get_gradients(Matrix& psi_grad, const int& n, 
    const std::vector<int>& row, const std::vector<int>& col) const;
private:
  std::unique_ptr<GroundState> groundstate_;
  //wf_descriptor wf_;
  std::string name_;
  //bool pairing_type_{false};
  //wf_type type_;
  unsigned num_sites_;
  // BCS_state bcs_state_;
  // FS_state fermisea_;
  Matrix psiup_;
  Matrix psidn_;

  std::vector<Matrix> psiup_grad_;
  std::vector<Matrix> psidn_grad_;
  bool single_determinant_{true};
  bool pairwf_{true};
  bool have_gradient_{false};

  // The great trick of accessing the wf matrix using pointers
#ifdef REAL_WAVEFUNCTION
  double* p_psiup_;
  double* p_psidn_;
  std::vector<double *> p_psiup_grad_;
  std::vector<double *> p_psidn_grad_;
#else
  std::complex<double>* p_psiup_;
  std::complex<double>* p_psidn_;
  std::vector<std::complex<double>*> p_psiup_grad_;
  std::vector<std::complex<double>*> p_psidn_grad_;
#endif
  // matrices & solvers
};


} // end namespace var

#endif
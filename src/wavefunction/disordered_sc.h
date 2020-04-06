/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-22 22:41:54
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 16:01:40
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef DISORDERED_SC_H
#define DISORDERED_SC_H

#include <Eigen/Sparse>
#include "./groundstate.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class MatrixElem
{
public:
  MatrixElem() : bond_phase_{1} {}
  MatrixElem(const unsigned& row, const unsigned& col, const double& val) 
    : row_{row}, col_{col}, value_{val}, bond_phase_{1} {}
  MatrixElem(const unsigned& row, const unsigned& col, const double& val, 
    const int& bphase) 
    : row_{row}, col_{col}, value_{val}, bond_phase_{bphase} {}
  ~MatrixElem() {}
  void change_value(const double& val) { value_=val; }
  void set_bond_phase(const int& bphase) { bond_phase_=bphase; }
  const unsigned& row(void) const { return row_; }
  const unsigned& col(void) const { return col_; }
  const double& value(void) const { return value_; }
  const int& bond_phase(void) const { return bond_phase_; }
private:
  unsigned row_{0};
  unsigned col_{0};
  double value_{0.0};
  int bond_phase_{1}; 
};

class DisorderedSC : public GroundState
{
public:
  DisorderedSC() : GroundState(true) {}
  DisorderedSC(const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph); 
  ~DisorderedSC() {} 
  int init(const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph);
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_gradient(std::vector<Matrix>& psi_gradient) override; 
private:
  double large_number_{1.0E+2};
  unsigned mu_start_;
  unsigned t_start_;
  unsigned delta_start_;
  // matrices

  using EigenTriplet = Eigen::Triplet<double>;
  std::vector<MatrixElem> quadratic_coeffs_; 
  std::vector<MatrixElem> pairing_coeffs_; 
  RealMatrix quadratic_ham_;
  RealMatrix pairing_ham_;
  RealMatrix work_;
  /*
  Eigen::SparseMatrix<double> quadratic_ham_;
  Eigen::SparseMatrix<double> pairing_ham_;
  Eigen::SparseMatrix<double> work_;
  ComplexMatrix psi_work2_;
  */
  RealMatrix delta_;
  std::vector<double> phi_k_;
  Matrix psi_work_;

  mutable Eigen::SelfAdjointEigenSolver<RealMatrix> es_quad_;
  mutable Eigen::SelfAdjointEigenSolver<RealMatrix> es_pair_;

  void hack_gradient(std::vector<Matrix>& psi_gradient); 
  void get_pair_amplitudes_sitebasis(const Eigen::VectorXd& es_eigenvalues, 
  const RealMatrix& es_eigenvectors, const Eigen::VectorXd& delta, Matrix& psi);
  void get_gradient_pairing_coeff(std::vector<Matrix>& psi_gradient,
    const unsigned& start_pos);
};


} // end namespace var


#endif
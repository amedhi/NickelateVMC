/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 14:51:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-22 00:35:51
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MF_MODEL_H
#define MF_MODEL_H

#include <iostream>
#include <string>
#include <complex>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../model/quantum_op.h"
#include "../model/model.h"
#include "./varparm.h"
#include "./blochbasis.h"
#include "./matrix.h"

enum class mf_order {none, af, fm, ssc, dsc, pplusip, af_dsc, disordered_sc};

//constexpr std::complex<double> ii(void) { return std::complex<double>{0.0,static_cast<double>(1.0)}; }

namespace var {
using vparm_t = std::pair<std::string,double>;

/*
//class SiteTerm : public std::unordered_map<unsigned, SiteOperatorTerm>
class Unitcell_Term
{
public:
  Unitcell_Term() {}
  ~Unitcell_Term() {}
  void build_bondterm(const model::HamiltonianTerm& sterm, const lattice::Lattice& lattice);
  void build_siteterm(const model::HamiltonianTerm& sterm, const lattice::Lattice& lattice);
  const unsigned& num_out_bonds(void) const { return num_out_bonds_; } 
  const Vector3d& bond_vector(const unsigned& i) const { return bond_vectors_[i]; }
  const ComplexMatrix& coeff_matrix(const unsigned& i=0) const { return coeff_matrices_[i]; }
  //const double& coupling(const unsigned& site_type) const; 
  const model::op::quantum_op& qn_operator(void) const { return op_; }
private:
  model::op::quantum_op op_;
  unsigned num_out_bonds_;
  std::vector<ComplexMatrix> coeff_matrices_;
  std::vector<Vector3d> bond_vectors_;
};*/

class UnitcellTerm
{
public:
  UnitcellTerm() {}
  ~UnitcellTerm() {}
  void build_bondterm(const model::HamiltonianTerm& sterm, const lattice::Lattice& lattice);
  void build_siteterm(const model::HamiltonianTerm& sterm, const lattice::Lattice& lattice);
  void eval_coupling_constant(const model::ModelParams& cvals, const model::ModelParams& pvals);
  const int& num_out_bonds(void) const { return num_out_bonds_; } 
  const Vector3d& bond_vector(const int& i) const { return bond_vectors_[i]; }
  const ComplexMatrix& coeff_matrix(const int& i=0) const { return coeff_matrices_[i]; }
  //const double& coupling(const unsigned& site_type) const; 
  const model::op::quantum_op& qn_operator(void) const { return op_; }
private:
  using strMatrix = std::vector<std::vector<std::string> >;
  model::op::quantum_op op_;
  int num_out_bonds_;
  int dim_;
  std::vector<ComplexMatrix> coeff_matrices_;
  std::vector<strMatrix> expr_matrices_;
  std::vector<Vector3d> bond_vectors_;
};

class MF_Model : public model::Hamiltonian
{
public:
  MF_Model() {}
  ~MF_Model() {}
  int init(const lattice::Lattice& lattice) override;
  int finalize(const lattice::Lattice& lattice);
  void update(const input::Parameters& inputs);
  void update_terms(void) override;
  void update_site_parameter(const std::string& pname, const double& pvalue);
  void construct_kspace_block(const Vector3d& kvec);
  const ComplexMatrix& quadratic_spinup_block(void) const { return quadratic_block_up_; }
  const ComplexMatrix& quadratic_spindn_block(void) const { return quadratic_block_dn_; }
  const ComplexMatrix& pairing_part(void) const { return pairing_block_; }

private:
  using Model = model::Hamiltonian;
  std::vector<UnitcellTerm> usite_terms_;
  std::vector<UnitcellTerm> ubond_terms_;
  // matrices in kspace representation
  int dim_;
  ComplexMatrix quadratic_block_up_;
  ComplexMatrix quadratic_block_dn_;
  ComplexMatrix pairing_block_;
  ComplexMatrix work, work2;

  void build_unitcell_terms(const lattice::Lattice& lattice);
  void update_unitcell_terms(void);
};


} // end namespace var

#endif
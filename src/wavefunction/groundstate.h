/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:32:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-13 10:47:01
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef GROUNDSTATE_H
#define GROUNDSTATE_H

#include <vector>
#include <Eigen/Eigenvalues>
#include "../scheduler/task.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

// Descriptor for measurement of SC correlation 
class CorrelationPairs
{
public:
  enum class corr_t {bond_singlet, site_singlet, null};
  CorrelationPairs() {}
 ~CorrelationPairs() {}
  void set_correlation_type(const corr_t& ctype) 
  { 
    clear_correlation_pairs();
    corr_type_ = ctype; 
  }
  void clear_correlation_pairs(void) { direction_bond_.clear(); corr_pairs_.clear(); }
  void add_correlation_sites(const int& direction_bond)
  {
    direction_bond_.push_back(direction_bond);
    corr_pairs_.push_back({-1,-1});
  }
  void add_correlation_bonds(const int& direction_bond, const std::pair<int,int>& bond_types)
  {
    direction_bond_.push_back(direction_bond);
    corr_pairs_.push_back(bond_types);
  }
  const corr_t& correlation_type(void) const { return corr_type_; }
  int num_correlation_pairs(void) const { return direction_bond_.size(); }
  const int& direction_bond(const int& i) const { return direction_bond_[i]; }
  const std::pair<int,int>& correlation_pair(const int& i) const { return corr_pairs_[i]; }
  const std::vector<std::pair<int,int>>& correlation_pairs(void) const { return corr_pairs_; }
private:
  corr_t corr_type_{corr_t::bond_singlet};
  std::vector<int> direction_bond_;
  std::vector<std::pair<int,int>> corr_pairs_;
};

class MF_Order : public CorrelationPairs 
{
public:
  enum class order_t {SC, AF, AFSC, CDW, null};
  enum class pairing_t {SWAVE, DWAVE, EXTENDED_S, DXY, D_PLUS_ID, CUSTOM, SC_CDW_SDW, null};
  MF_Order() {}
  MF_Order(const MF_Order::order_t& order, const MF_Order::pairing_t& pairsymm) 
    : order_{order}, pairing_{pairsymm} {}
  ~MF_Order() {} 
  int id(void) const;
  void set_order(const MF_Order::order_t& order) { order_=order; }
  void set_pair_symmetry(const MF_Order::pairing_t& pairsymm) { pairing_=pairsymm; }
  const MF_Order::order_t& order(void) const { return order_; }
  const MF_Order::pairing_t& pair_symm(void) const { return pairing_; } 
  bool pairing_type(void) const 
  { 
    if (pairing_==MF_Order::pairing_t::null) return false;
    else return true;
  }
  //std::vector<std::pair<int,int>>& correlation_pairs(void) { return correlation_pairs_; }
  //const std::vector<std::pair<int,int>>& correlation_pairs(void) const { return correlation_pairs_; }
private:
  order_t order_{order_t::null};
  pairing_t pairing_{pairing_t::null};
  std::vector<std::pair<int,int>> correlation_pairs_;
};

class GroundState : public MF_Order
{
public:
  GroundState() : MF_Order() {}
  GroundState(const MF_Order::order_t& order, const MF_Order::pairing_t& pair_symm)
    : MF_Order(order, pair_symm) {}
  virtual ~GroundState() {} 
  virtual void update(const input::Parameters& inputs);
  virtual void update(const var::parm_vector& pvector, const unsigned& start_pos=0);
  virtual void update(const lattice::Lattice& lattice);
  virtual void get_wf_amplitudes(Matrix& psi);
  virtual void get_wf_amplitudes(Matrix& psiup, Matrix& psidn);
  virtual void get_wf_gradient(std::vector<Matrix>& psi_grad); 
  virtual void get_wf_gradient(std::vector<Matrix>& psiup_grad, std::vector<Matrix>& psidn_grad); 
  virtual std::string info_str(void) const; 
  const VariationalParms& varparms(void) { return varparms_; }
  const int& num_varparms(void) const { return num_varparms_; }
  const bool& is_nonmagnetic(void) const { return nonmagnetic_; }
  const int& num_sites(void) const { return num_sites_; }
  const int& num_spins(void) const { return num_spins_; }
  const int& num_upspins(void) const { return num_upspins_; }
  const int& num_dnspins(void) const { return num_dnspins_; }
  const double& hole_doping(void) const { return hole_doping_; }
  const basis::BlochBasis blochbasis(void) const { return blochbasis_; }
protected:
  std::string name_;
  int num_sites_{0};
  int num_bonds_{0};
  int num_kpoints_{0};
  int kblock_dim_{0};
  int num_varparms_{0};
  basis::BlochBasis blochbasis_;
  MF_Model mf_model_;
  MF_Order mf_order_;
  VariationalParms varparms_;
  // fourier transform matrix
  ComplexMatrix FTU_;
  // solvers
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_k_up;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_minusk_dn;
  void set_nonmagnetic(const bool& yesno) { nonmagnetic_=yesno; }
  void set_particle_num(const input::Parameters& inputs);
  void reset_spin_num(const int& num_upspin, const int& num_dnspin);
  void set_ft_matrix(const lattice::Lattice& lattice);
  double get_noninteracting_mu(void);
private:
  bool nonmagnetic_{true};
  int num_spins_{0};
  int num_upspins_{0};
  int num_dnspins_{0};
  double last_hole_doping_{2.5}; // unlikely input
  double hole_doping_inp_{0.0};
  double hole_doping_{0.0};
  double band_filling_{1.0};
  //MF_Model mf_model_;
  //basis::BlochBasis blochbasis_;
};


} // end namespace var

#endif
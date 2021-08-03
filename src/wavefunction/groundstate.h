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

class MF_Order 
{
public:
  enum class order_t {SC, AF, AFSC, CDW, null};
  enum class pairing_t {SWAVE, DWAVE, EXTENDED_S, D_PLUS_ID, CUSTOM, null};
  MF_Order() {}
  MF_Order(const MF_Order::order_t& order, const MF_Order::pairing_t& pairsymm) 
    : order_{order}, pairing_{pairsymm} {}
  ~MF_Order() {} 
  int id(void) const;
  void set_order(const MF_Order::order_t& order) { order_=order; }
  void set_pair_symmetry(const MF_Order::pairing_t& pairsymm) { pairing_=pairsymm; }
  const MF_Order::order_t& order(void) const { return order_; }
  const MF_Order::pairing_t& pair_symm(void) { return pairing_; } 
  bool pairing_type(void) const 
  { 
    if (pairing_==MF_Order::pairing_t::null) return false;
    else return true;
  }
private:
  order_t order_{order_t::null};
  pairing_t pairing_{pairing_t::null};
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
  virtual void update(const lattice::LatticeGraph& graph);
  virtual void get_wf_amplitudes(Matrix& psi);
  virtual void get_wf_gradient(std::vector<Matrix>& psi_gradient); 
  virtual std::string info_str(void) const; 
  const VariationalParms& varparms(void) { return varparms_; }
  const int& num_varparms(void) const { return num_varparms_; }
  const bool& if_nonmagnetic(void) const { return nonmagnetic_; }
  const int& num_upspins(void) const { return num_upspins_; }
  const int& num_dnspins(void) const { return num_dnspins_; }
  const double& hole_doping(void) const { return hole_doping_; }
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
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_minusk_up;
  void set_nonmagnetic(const bool& yesno) { nonmagnetic_=yesno; }
  void set_particle_num(const input::Parameters& inputs);
  void set_ft_matrix(const lattice::LatticeGraph& graph);
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
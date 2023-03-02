/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 22:41:38
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-13 15:03:27
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef BCS_STATE_H
#define BCS_STATE_H

#include "./groundstate.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class BCS_State : public GroundState
{
public:
  BCS_State() {}
  BCS_State(const MF_Order::order_t& order, const MF_Order::pairing_t& pair_symm,
    const input::Parameters& inputs, const lattice::Lattice& lattice,
    const model::Hamiltonian& model); 
  virtual ~BCS_State() {} 
  int init(const input::Parameters& inputs, 
    const lattice::Lattice& lattice, const model::Hamiltonian& model); 
  std::string info_str(void) const override; 
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void update(const lattice::Lattice& lattice) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_gradient(std::vector<Matrix>& psi_gradient) override; 
private:
  std::string order_name_;
  bool interband_pairing_{false};
  bool wf_analytical_form_{false};
  bool noninteracting_mu_{true};
  double singular_ampl_{1.0E+2};
  lattice::lattice_id lattice_id_;
  // matrices
  ComplexMatrix work_;
  ComplexMatrix delta_k_;
  ComplexMatrix dphi_k_;
  ComplexMatrix bdg_mat_;
  ComplexMatrix uk_;
  ComplexMatrix vk_;
  ComplexMatrix ubk_;
  ComplexMatrix vbk_;
  std::vector<ComplexMatrix> phi_k_;
  std::vector<ComplexMatrix> work_k_;
  mutable Eigen::SelfAdjointEigenSolver<ComplexMatrix> es_bdg_;

  //void add_chemical_potential(const input::Parameters& inputs);
  void get_pair_amplitudes_oneband(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_intraband(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_interband(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, Matrix& psi);
  void analytical_amplitudes_NICKELATE2L(std::vector<ComplexMatrix>& phi_k);
  void analytical_gradient_NICKELATE2L(std::vector<Matrix>& psi_gradient);
  double get_mf_energy(void);
};


} // end namespace var


#endif
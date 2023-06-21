/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-02-20 12:21:42
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-02-20 12:21:42
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef FERMISEA_H
#define FERMISEA_H

#include "./groundstate.h"
#include "./mf_model.h"
#include "./matrix.h"

namespace var {

class Fermisea : public GroundState
{
public:
  Fermisea() : GroundState() {}
  Fermisea(const MF_Order::order_t& order, const input::Parameters& inputs, 
    const lattice::Lattice& lattice, const model::Hamiltonian& model);
  ~Fermisea() {} 
  int init(const input::Parameters& inputs, const lattice::Lattice& lattice,
    const model::Hamiltonian& model);
  std::string info_str(void) const override; 
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void update(const lattice::Lattice& lattice) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_amplitudes(Matrix& psiup, Matrix& psidn) override;
  void get_wf_gradient(std::vector<Matrix>& psi_grad) override; 
  void get_wf_gradient(std::vector<Matrix>& psiup_grad, std::vector<Matrix>& psidn_grad) override; 
private:
  std::string order_name_{"NULL"};
  bool mu_variational_{false};
  bool degeneracy_warning_{false};
  lattice::lattice_id lattice_id_;
  // ground state
  //bool have_TP_symmetry_{true};
  double fermi_energy_;
  double total_energy_;
  struct kshell_t {int k; int nmin; int nmax;};
  std::vector<kshell_t> kshells_up_;
  std::vector<kshell_t> kshells_dn_;

  // SC correlaton function
  int rmax_;
  Vector3d alpha_;
  Vector3d beta_;
  std::vector<Vector3d> R_list_;
  Eigen::MatrixXcd corr_aa_;
  Eigen::MatrixXcd corr_ab_;
  Eigen::MatrixXcd corr_fs_;

  void construct_groundstate(void);
  void get_amplitudes_sitebasis(ComplexMatrix& psiup, ComplexMatrix& psidn);
  void get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k);
  void get_sc_correlation(void);
  double get_mf_energy(void);
};


} // end namespace var


#endif

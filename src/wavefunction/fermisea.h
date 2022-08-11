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
    const lattice::Lattice& lattice); 
  ~Fermisea() {} 
  int init(const input::Parameters& inputs, const lattice::Lattice& lattice);
  std::string info_str(void) const override; 
  void update(const input::Parameters& inputs) override;
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0) override;
  void update(const lattice::Lattice& lattice) override;
  void get_wf_amplitudes(Matrix& psi) override;
  void get_wf_gradient(std::vector<Matrix>& psi_gradient) override; 
private:
  std::string order_name_{"NULL"};
  bool noninteracting_mu_{true};
  // ground state
  bool have_TP_symmetry_{true};
  double fermi_energy_;
  double total_energy_;
  bool degeneracy_warning_{false};
  struct kshell_t {int k; int nmin; int nmax;};
  std::vector<kshell_t> kshells_up_;
  std::vector<kshell_t> kshells_dn_;

  // matrices
  ComplexMatrix work_;
  std::vector<ComplexMatrix> phi_k_;

  // SC correlaton function
  int rmax_;
  Vector3d alpha_;
  Vector3d beta_;
  std::vector<Vector3d> R_list_;
  Eigen::MatrixXcd corr_aa_;
  Eigen::MatrixXcd corr_ab_;
  Eigen::MatrixXcd corr_fs_;

  void add_chemical_potential(const input::Parameters& inputs);
  void construct_groundstate(void);
  void get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k);
  void get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, Matrix& psi);
  void get_sc_correlation(void);
  double get_mf_energy(void);
};


} // end namespace var


#endif

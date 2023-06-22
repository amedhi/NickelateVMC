/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-18 13:54:54
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-10 23:35:53
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef SYSCONFIG_H
#define SYSCONFIG_H

#include <memory>
#include "../scheduler/worker.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "../wavefunction/wavefunction.h"
#include "../wavefunction/projector.h"
#include "../wavefunction/matrix.h"
#include "./basisstate.h"

namespace vmc {

constexpr double dratio_cutoff(void) { return 1.0E-8; } 
constexpr double gfactor_cutoff(void) { return 1.0E-8; } 

class SysConfig : public BasisState
{
public:
  SysConfig(const input::Parameters& parms, const lattice::Lattice& lattice, 
    const model::Hamiltonian& model);
  ~SysConfig() {}
  std::string info_str(void) const; 
  int build(const lattice::Lattice& lattice, const input::Parameters& inputs, 
    const bool& with_gradient=false);
  int build(const lattice::Lattice& lattice, const var::parm_vector& vparms, 
    const bool& need_psi_grad=false);
  int rebuild(const lattice::Lattice& lattice);
  std::string signature_str(void) const { return wf.signature_str(); } 
  const unsigned& num_varparms(void) const { return num_varparms_; } 
  int num_particles(void) const { return num_upspins_+num_dnspins_; } 
  const var::parm_vector& vparm_values(void);
  const std::vector<double>& vparm_vector(void); 
  const std::vector<std::string>& varp_names(void) const { return vparm_names_; }
  const var::parm_vector& vparm_lbound(void) const { return vparm_lbound_; } 
  const var::parm_vector& vparm_ubound(void) const { return vparm_ubound_; } 
  const double& hole_doping(void) const { return wf.hole_doping(); }
  int update_state(void);
  double accept_ratio(void);
  void reset_accept_ratio(void);
  int apply_ni_dblon(const int& site_i) const;
  int apply_ni_holon(const int& site_i) const;
  int apply(const model::op::quantum_op& qn_op, const int& site_i) const;
  amplitude_t apply(const model::op::quantum_op& op, const int& fr_site, 
    const int& to_site, const int& bc_state, const std::complex<double>& bc_phase) const;
  amplitude_t apply_bondsinglet_hop(const int& fr_site_i, 
    const int& fr_site_ia, const int& to_site_j, const int& to_site_jb) const;
  amplitude_t apply_sitepair_hop(const int& fr_site, const int& to_site) const;
  amplitude_t apply_upspin_hop(const int& i, const int& j,
    const int& bc_state, const std::complex<double>& bc_phase) const;
  amplitude_t apply_dnspin_hop(const int& i, const int& j,
    const int& bc_state, const std::complex<double>& bc_phase) const;
  amplitude_t apply_cdagc_up(const int& fr_site, const int& to_site,
    const int& bc_state, const std::complex<double>& bc_phase) const;
  amplitude_t apply_cdagc_dn(const int& fr_site, const int& to_site,
    const int& bc_state, const std::complex<double>& bc_phase) const;
  amplitude_t apply_sisj_plus(const int& i, const int& j) const;
  void get_grad_logpsi(RealVector& grad_logpsi) const;
  const int& num_updates(void) const { return num_updates_; }
  const var::Wavefunction& wavefunc(void) const { return wf; }
  void print_stats(std::ostream& os=std::cout) const;
  //var::VariationalParms& var_parms(void) { return wf.var_parms(); }
private:
  bool single_determinant_{true};
  var::Wavefunction wf;
  var::WavefunProjector pj;
  Matrix psiup_mat_;
  Matrix psidn_mat_;
  Matrix psiup_inv_;
  Matrix psidn_inv_;
  mutable Matrix psiup_grad_;
  mutable Matrix psidn_grad_;

  mutable ColVector psi_row;
  mutable RowVector psi_col;
  mutable RowVector inv_row;
  int num_sites_;
  int num_upspins_;
  int num_dnspins_;

  // variational parameters
  unsigned num_pj_parms_{0};
  unsigned num_wf_parms_{0};
  unsigned num_varparms_{0};
  mutable var::parm_vector vparm_values_;
  mutable std::vector<double> vparm_vector_;
  std::vector<std::string> vparm_names_;
  var::parm_vector vparm_lbound_;
  var::parm_vector vparm_ubound_;

  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_updates_{0};
  //int num_total_steps_{0};
  int num_uphop_moves_{0};
  int num_dnhop_moves_{0};
  int num_exchange_moves_{0};
  int num_iterations_{0};
  int refresh_cycle_{100};
  long num_proposed_moves_[move_t::end];
  long num_accepted_moves_[move_t::end];
  long last_proposed_moves_;
  long last_accepted_moves_;

  // helper methods
  int init_config(void);
  int set_run_parameters(void);
  int do_upspin_hop(void);
  int do_dnspin_hop(void);
  int do_spin_exchange(void);
  int inv_update_upspin(const int& upspin, const ColVector& psi_row, 
    const amplitude_t& det_ratio);
  int inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
    const amplitude_t& det_ratio);
};

} // end namespace vmc

#endif
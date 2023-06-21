/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-18 14:01:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:16:31
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "./sysconfig.h"
#include <Eigen/SVD>

namespace vmc {

SysConfig::SysConfig(const input::Parameters& inputs, 
  const lattice::Lattice& lattice, const model::Hamiltonian& model)
  : BasisState(lattice, model)
  , wf(lattice, inputs, model)
  , pj(lattice, inputs)
  , num_sites_(lattice.num_sites())
{
  // set hubbard sites 
  /*
  if (lattice.id()==lattice::lattice_id::NICKELATE ||
      lattice.id()==lattice::lattice_id::NICKELATE_2D ||
      lattice.id()==lattice::lattice_id::NICKELATE_2L) {
    if (model.id()==model::model_id::HUBBARD) {
      for (const auto& s : lattice.sites()) {
        if (s.type() == 1) set_projection(s.id(), proj_t::null);
        else set_projection(s.id(), proj_t::partial);
      }
    }
    else if (model.id()==model::model_id::TJ) {
      for (const auto& s : lattice.sites()) {
        if (s.type() == 1) set_projection(s.id(), proj_t::null);
        else set_projection(s.id(), proj_t::full);
      }
    }
    else {
      throw std::range_error("*error: sysconfig: internal error");
    }
  }
  */

  // If model has full projection, then switch off GW projector
  if (model.projection()) {
    pj.gw_switch_off();
  }

  // particle number
  num_upspins_ = wf.num_upspins();
  num_dnspins_ = wf.num_dnspins();

  // variational parameters
  num_pj_parms_ = pj.varparms().size();
  num_wf_parms_ = wf.varparms().size();
  num_varparms_ = (num_pj_parms_ + num_wf_parms_);
  vparm_names_.resize(num_varparms_);
  vparm_lbound_.resize(num_varparms_);
  vparm_ubound_.resize(num_varparms_);
  // names
  pj.get_vparm_names(vparm_names_,0);
  wf.get_vparm_names(vparm_names_,num_pj_parms_);
  // values are not static and may change
  // bounds
  pj.get_vparm_lbound(vparm_lbound_,0);
  wf.get_vparm_lbound(vparm_lbound_,num_pj_parms_);
  pj.get_vparm_ubound(vparm_ubound_,0);
  wf.get_vparm_ubound(vparm_ubound_,num_pj_parms_);
}

const var::parm_vector& SysConfig::vparm_values(void) 
{
  // values as 'var::parm_vector'
  vparm_values_.resize(num_varparms_);
  pj.get_vparm_values(vparm_values_,0);
  wf.get_vparm_values(vparm_values_,num_pj_parms_);
  return vparm_values_;
}

const std::vector<double>& SysConfig::vparm_vector(void) 
{
  // values as 'std::double'
  vparm_vector_.resize(num_varparms_);
  pj.get_vparm_vector(vparm_vector_,0);
  wf.get_vparm_vector(vparm_vector_,num_pj_parms_);
  return vparm_vector_;
}

std::string SysConfig::info_str(void) const
{
  std::ostringstream info;
  info << wf.info_str();
  info.precision(6);
  info.setf(std::ios_base::fixed);
  info << "# Variational parameters:\n";
  for (const auto& v : wf.varparms()) 
    info << "# " << v.name() << " = " << v.value() << "\n";
  for (const auto& v : pj.varparms()) 
    info << "# " << v.name() << " = " << v.value() << "\n";
  return info.str();
}

int SysConfig::build(const lattice::Lattice& lattice, const input::Parameters& inputs,
    const bool& with_gradient)
{
  if (num_sites_==0) return -1;
  pj.update(inputs);
  wf.compute(lattice, inputs, with_gradient);
  return init_config();
}

int SysConfig::build(const lattice::Lattice& lattice, const var::parm_vector& pvector,
  const bool& need_psi_grad)
{
  if (num_sites_==0) return -1;
  pj.update(pvector, 0);
  unsigned start_pos = pj.varparms().size();
  wf.compute(lattice, pvector, start_pos, need_psi_grad);
  return init_config();
}

// rebuild for new lattice boundary twist
int SysConfig::rebuild(const lattice::Lattice& lattice)
{
  wf.recompute(lattice);
  return init_config();
}


int SysConfig::init_config(void)
{
  num_upspins_ = wf.num_upspins();
  num_dnspins_ = wf.num_dnspins();
  if (num_upspins_==0 && num_dnspins_==0) return -1;
  if (num_upspins_ != num_dnspins_) 
    throw std::range_error("*SysConfig::init_config: unequal UP & DN spin case not implemented");
  /*
  // small 'gfactor' caution
  bool tmp_restriction = false;
  bool original_state = BasisState::double_occupancy();
  if (pj.gw_projection()) {
    if (pj.gw_projection_strong()) {
      BasisState::allow_double_occupancy(false);
      tmp_restriction = true;
    }
  }*/

  BasisState::init_spins(num_upspins_, num_dnspins_);
  psiup_mat_.resize(num_upspins_, num_dnspins_);
  psi_inv.resize(num_dnspins_, num_upspins_);
  // try for a well condictioned amplitude matrix
  double rcond = 0.0;
  int num_attempt = 0;
  while (rcond<1.0E-20) {
  //while (rcond<1.0E-12) {
    BasisState::set_random();
    //BasisState::set_custom();
    wf.get_amplitudes(psiup_mat_,upspin_sites(),dnspin_sites());
    //std::cout << "psiup_mat_ = \n";
    //std::cout << psiup_mat_ << "\n"; getchar();
    // reciprocal conditioning number
    Eigen::JacobiSVD<Matrix> svd(psiup_mat_);
    // reciprocal cond. num = smallest eigenval/largest eigen val
    rcond = svd.singularValues()(svd.singularValues().size()-1)/svd.singularValues()(0);
    if (std::isnan(rcond)) rcond = 0.0; 
    //std::cout << "rcondition number = "<< rcond << "\n"; getchar();
    if (++num_attempt > 1000) {
      throw std::underflow_error("*SysConfig::init: configuration wave function ill conditioned.");
    }
  }
  //if (tmp_restriction) allow_double_occupancy(original_state);

  // amplitude matrix invers
  psi_inv = psiup_mat_.inverse();
  // check
  //std::cout << psiup_mat_ << "\n"; getchar();
  //std::cout << psi_inv << "\n"; getchar();

  // run parameters
  set_run_parameters();
  return 0;
}

int SysConfig::set_run_parameters(void)
{
  num_updates_ = 0;
  num_iterations_ = 0;
  refresh_cycle_ = 100;
  // number of moves per mcstep
  int n_up = static_cast<int>(num_upspins_);
  int n_dn = static_cast<int>(num_dnspins_);
  //if (double_occupancy()) {
  if (true) {
    num_uphop_moves_ = num_upspins_;
    num_dnhop_moves_ = num_dnspins_;
    num_exchange_moves_ = std::min(n_up, n_dn);
    //num_exchange_moves_ = 2*std::min(n_up, n_dn);
  }
  else {
    int num_holes = num_sites_-(num_upspins_+num_dnspins_);
    num_uphop_moves_ = std::min(n_up, num_holes);
    num_dnhop_moves_ = std::min(n_dn, num_holes);
    num_exchange_moves_ = std::min(n_up, n_dn);
    //num_exchange_moves_ = 4*std::min(n_up, n_dn);
  }
  for (int i=0; i<move_t::end; ++i) {
    num_proposed_moves_[i] = 0;
    num_accepted_moves_[i] = 0;
  }
  last_proposed_moves_ = 1;
  last_accepted_moves_ = 1;

  // work arrays 
  psi_row.resize(num_dnspins_);
  psi_col.resize(num_upspins_);
  inv_row.resize(num_upspins_);
  psi_grad.resize(num_upspins_,num_dnspins_);
  return 0;
}

int SysConfig::update_state(void)
{
  for (int n=0; n<num_uphop_moves_; ++n) do_upspin_hop();
  for (int n=0; n<num_dnhop_moves_; ++n) do_dnspin_hop();
  for (int n=0; n<num_exchange_moves_; ++n) do_spin_exchange();
  num_updates_++;
  num_iterations_++;

  //auto psi = psiup_mat_.determinant();
  //std::cout<<psi.real()<<"   "<<psi.imag()<<"\n";
  //std::cout << psi_inv << "\n";

  if (num_iterations_ == refresh_cycle_) {
    psi_inv = psiup_mat_.inverse();
    num_iterations_ = 0;
  }
  return 0;
}

//int SysConfig::dh_ratio()

int SysConfig::do_upspin_hop(void)
{
  if (!gen_upspin_hop()) return 1;
  num_proposed_moves_[move_t::uphop]++;
  last_proposed_moves_++;
  int upspin = which_upspin();
  int to_site = which_site();
  // new row for this move
  wf.get_amplitudes(psi_row, to_site, dnspin_sites());
  amplitude_t det_ratio = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  //std::cout << "psi_inv = " << psi_inv.col(upspin).transpose() << "\n"; 
  //std::cout << "det_ratio = " << det_ratio << "\n"; //getchar();
  if (std::abs(det_ratio) < dratio_cutoff()) return 0; // for safety
  //double proj_ratio = pj.gw_ratio(dblocc_increament());
  double proj_ratio = pj.gw_ratio(*this, which_frsite(), to_site);

  amplitude_t weight_ratio = det_ratio * proj_ratio;
  double transition_proby = std::norm(weight_ratio);
  //std::cout << "W = " << transition_proby << "\n"; getchar();
  if (rng().random_real()<transition_proby) {
    num_accepted_moves_[move_t::uphop]++;
    last_accepted_moves_++;
    // upddate state
    accept_last_move();
    // update amplitudes
    inv_update_upspin(upspin,psi_row,det_ratio);
  }
  return 0;
}

int SysConfig::do_dnspin_hop(void)
{
  if (!gen_dnspin_hop()) return 1;
  num_proposed_moves_[move_t::dnhop]++;
  last_proposed_moves_++;
  int dnspin = which_dnspin();
  int to_site = which_site();
  // new col for this move
  wf.get_amplitudes(psi_col, upspin_sites(), to_site);
  amplitude_t det_ratio = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
  if (std::abs(det_ratio) < dratio_cutoff()) return 0; // for safety
  //double proj_ratio = pj.gw_ratio(dblocc_increament());
  double proj_ratio = pj.gw_ratio(*this, which_frsite(), to_site);

  amplitude_t weight_ratio = det_ratio * proj_ratio;
  double transition_proby = std::norm(weight_ratio);
  if (rng().random_real()<transition_proby) {
    num_accepted_moves_[move_t::dnhop]++;
    last_accepted_moves_++;
    // upddate state
    accept_last_move();
    // update amplitudes
    inv_update_dnspin(dnspin,psi_col,det_ratio);
  }
  return 0;
}

int SysConfig::do_spin_exchange(void)
{
  if (!gen_exchange_move()) return 1;
  int upspin = which_upspin();
  int up_tosite = which_upsite();
  num_proposed_moves_[move_t::exch]++;
  last_proposed_moves_++;
  // for upspin hop forward
  wf.get_amplitudes(psi_row, up_tosite, dnspin_sites());
  amplitude_t det_ratio1 = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  if (std::abs(det_ratio1) < dratio_cutoff()) return 0; // for safety

  // now for dnspin hop backward
  int dnspin = which_dnspin();
  int dn_tosite = which_dnsite();
  // new col for this move
  wf.get_amplitudes(psi_col, upspin_sites(), dn_tosite);
  // since the upspin should have moved
  wf.get_amplitudes(psi_col(upspin), up_tosite, dn_tosite);
  // updated 'dnspin'-th row of psi_inv
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;
  // elements other than 'upspin'-th
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
  }
  for (int i=upspin+1; i<static_cast<int>(num_upspins_); ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
  }
  inv_row(upspin) = ratio_inv * psi_inv(dnspin,upspin);
  // ratio for the dnspin hop
  amplitude_t det_ratio2 = psi_col.cwiseProduct(inv_row).sum();
  if (std::abs(det_ratio2) < dratio_cutoff()) return 0; // for safety

  // weight ratio (gw pj does not play here)
  amplitude_t weight_ratio = det_ratio1 * det_ratio2;
  double transition_proby = std::norm(weight_ratio);
  //std::cout << "W = " << transition_proby << "\n";
  if (rng().random_real()<transition_proby) {
    num_accepted_moves_[move_t::exch]++;
    last_accepted_moves_++;
    // upddate state
    accept_last_move();
    // update amplitudes
    inv_update_upspin(upspin,psi_row,det_ratio1);
    inv_update_dnspin(dnspin,psi_col,det_ratio2);
  }
  return 0;
}

int SysConfig::inv_update_upspin(const int& upspin, const ColVector& psi_row, 
  const amplitude_t& det_ratio)
{
  psiup_mat_.row(upspin) = psi_row;
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    psi_inv.col(i) -= beta * psi_inv.col(upspin);
  }
  for (int i=upspin+1; i<num_upspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    psi_inv.col(i) -= beta * psi_inv.col(upspin);
  }
  psi_inv.col(upspin) *= ratio_inv;
  return 0;
}

int SysConfig::inv_update_dnspin(const int& dnspin, const RowVector& psi_col, 
  const amplitude_t& det_ratio)
{
  psiup_mat_.col(dnspin) = psi_col;
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio;
  for (int i=0; i<dnspin; ++i) {
    amplitude_t beta = ratio_inv*psi_col.cwiseProduct(psi_inv.row(i)).sum();
    psi_inv.row(i) -= beta * psi_inv.row(dnspin);
  }
  for (int i=dnspin+1; i<num_dnspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_col.cwiseProduct(psi_inv.row(i)).sum();
    psi_inv.row(i) -= beta * psi_inv.row(dnspin);
  }
  psi_inv.row(dnspin) *= ratio_inv;
  return 0;
}

int SysConfig::apply_ni_dblon(const int& site_i) const
{
  return op_ni_dblon(site_i);
}

int SysConfig::apply_ni_holon(const int& site_i) const
{
  return op_ni_holon(site_i);
}

int SysConfig::apply(const model::op::quantum_op& qn_op, const int& site_i) const
{
  switch (qn_op.id()) {
    case model::op_id::ni_sigma:
      return op_ni_up(site_i)+op_ni_dn(site_i);
    case model::op_id::ni_up:
      return op_ni_up(site_i);
    case model::op_id::ni_dn:
      return op_ni_dn(site_i);
    case model::op_id::niup_nidn:
      return op_ni_updn(site_i); 
    case model::op_id::Sz:
      return op_Sz(site_i); 
    default: 
      throw std::range_error("SysConfig::apply: undefined site operator");
  }
}

amplitude_t SysConfig::apply(const model::op::quantum_op& qn_op, const int& fr_site, 
    const int& to_site, const int& bc_state, const std::complex<double>& bc_phase) const
{
  amplitude_t term(0); 
  switch (qn_op.id()) {
    case model::op_id::cdagc_sigma:
      term = apply_upspin_hop(fr_site,to_site,bc_state,bc_phase);
      term+= apply_dnspin_hop(fr_site,to_site,bc_state,bc_phase);
      break;
    case model::op_id::sisj_plus:
      term = apply_sisj_plus(fr_site,to_site); 
      break;
    case model::op_id::ni_nj:
      term = op_ni(fr_site)*op_ni(to_site); 
      break;
    default: 
      throw std::range_error("SysConfig::apply: undefined bond operator.");
  }
  return term;
}

amplitude_t SysConfig::apply_upspin_hop(const int& site_i, const int& site_j,
  const int& bc_state, const std::complex<double>& bc_phase) const
{
  if (site_i == site_j) return ampl_part(op_ni_up(site_i));
  if (op_cdagc_up(site_i,site_j)) {
    int upspin = which_upspin();
    //int delta_nd = dblocc_increament();
    int to_site = site_j;
    // det_ratio for the term
    wf.get_amplitudes(psi_row, to_site, dnspin_sites());
    amplitude_t det_ratio = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
    if (std::isnan(det_ratio)) {
      std::cout << site_i << "\n";
      std::cout << site_j << "\n";
      //std::cout << psi_row.transpose() << "\n";
      std::cout << psi_inv.col(upspin).transpose() << "\n";
      std::cout << det_ratio << "\n";
      getchar();
    }
    if (bc_state < 0) { // it's a boundary bond
      //return ampl_part(bc_phase*std::conj(det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(bc_phase*std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
    else {
      //return ampl_part(std::conj(det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
  }
  else if (op_cdagc_up(site_j,site_i)) {
    int upspin = which_upspin();
    //int delta_nd = dblocc_increament();
    int to_site = site_i;
    // det_ratio for the term
    wf.get_amplitudes(psi_row, to_site, dnspin_sites());
    amplitude_t det_ratio = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
    if (bc_state < 0) { // it's a boundary bond
      //return ampl_part(std::conj(bc_phase*det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(bc_phase*det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
    else {
      //return ampl_part(std::conj(det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
  }
  else {
    return amplitude_t(0.0);
  }

  /*
  if (i == j) return ampl_part(op_ni_up(i));
  if (!op_cdagc_up_plus(i,j)) return amplitude_t(0.0);
  int upspin = which_upspin();
  int to_site = which_site();
  int delta_nd = dblocc_increament();
  // det_ratio for the term
  wf.get_amplitudes(psi_row, to_site, dnspin_sites());
  amplitude_t det_ratio = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  det_ratio = ampl_part(std::conj(det_ratio));
  return amplitude_t(bc_phase) * det_ratio * pj.gw_ratio(delta_nd);
  */

}

amplitude_t SysConfig::apply_cdagc_up(const int& fr_site, const int& to_site,
  const int& bc_state, const std::complex<double>& bc_phase) const
{
  if (fr_site == to_site) return ampl_part(op_ni_up(fr_site));
  if (!op_cdagc_up(fr_site,to_site)) return amplitude_t(0.0);
  int upspin = which_upspin();
  //int delta_nd = dblocc_increament();
  // det_ratio for the term
  wf.get_amplitudes(psi_row, to_site, dnspin_sites());
  amplitude_t det_ratio = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  //det_ratio = ampl_part(std::conj(det_ratio));
  if (bc_state < 0) { // it's a boundary bond
    //return ampl_part(bc_phase*std::conj(det_ratio))*pj.gw_ratio(delta_nd);
    return ampl_part(bc_phase*std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
  }
  else {
    //return ampl_part(std::conj(det_ratio))*pj.gw_ratio(delta_nd);
    return ampl_part(std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
  }
}

amplitude_t SysConfig::apply_dnspin_hop(const int& site_i, const int& site_j,
  const int& bc_state, const std::complex<double>& bc_phase) const
{
  if (site_i == site_j) return ampl_part(op_ni_dn(site_i));
  if (op_cdagc_dn(site_i,site_j)) {
    int dnspin = which_dnspin();
    int to_site = site_j;
    //int delta_nd = dblocc_increament();
    // det_ratio for the term
    wf.get_amplitudes(psi_col, upspin_sites(), to_site);
    amplitude_t det_ratio = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
    if (bc_state < 0) { // it's a boundary bond
      //return ampl_part(std::conj(bc_phase*det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(bc_phase*det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
    else {
      //return ampl_part(std::conj(det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
  }
  else if (op_cdagc_dn(site_j,site_i)) {
    int dnspin = which_dnspin();
    int to_site = site_i;
    //int delta_nd = dblocc_increament();
    // det_ratio for the term
    wf.get_amplitudes(psi_col, upspin_sites(), to_site);
    amplitude_t det_ratio = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
    if (bc_state < 0) { // it's a boundary bond
      //return ampl_part(std::conj(bc_phase*det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(bc_phase*det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
    else {
      //return ampl_part(std::conj(det_ratio))*pj.gw_ratio(delta_nd);
      return ampl_part(std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
    }
  }
  else {
    return amplitude_t(0.0);
  }

  /*
  if (i == j) return ampl_part(op_ni_dn(i));
  if (!op_cdagc_dn_plus(i,j)) return amplitude_t(0.0);
  int dnspin = which_dnspin();
  int to_site = which_site();
  int delta_nd = dblocc_increament();
  // det_ratio for the term
  wf.get_amplitudes(psi_col, upspin_sites(), to_site);
  amplitude_t det_ratio = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
  det_ratio = ampl_part(std::conj(det_ratio));
  return amplitude_t(bc_phase) * det_ratio * pj.gw_ratio(delta_nd);
  */
}

amplitude_t SysConfig::apply_cdagc_dn(const int& fr_site, const int& to_site,
  const int& bc_state, const std::complex<double>& bc_phase) const
{
  if (fr_site == to_site) return ampl_part(op_ni_dn(fr_site));
  if (!op_cdagc_dn(fr_site,to_site)) return amplitude_t(0.0);
  int dnspin = which_dnspin();
  //int delta_nd = dblocc_increament();
  // det_ratio for the term
  wf.get_amplitudes(psi_col, upspin_sites(), to_site);
  amplitude_t det_ratio = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
  //det_ratio = ampl_part(std::conj(det_ratio));
  if (bc_state < 0) { // it's a boundary bond
    //return ampl_part(bc_phase*std::conj(det_ratio))*pj.gw_ratio(delta_nd);
    return ampl_part(bc_phase*std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
  }
  else {
    //return ampl_part(std::conj(det_ratio))*pj.gw_ratio(delta_nd);
    return ampl_part(std::conj(det_ratio))*pj.gw_ratio(*this,which_frsite(),to_site);
  }
}

amplitude_t SysConfig::apply_sisj_plus(const int& i, const int& j) const
{
/*----------------------------------------------------------------------
* Evaluates the following operator:
*   J = (S_i.S_j - (n_i n_j)/4)
* The operator can be cast in the form,
*   J = P_{ud} + P_{du}
* where,
*   P_{ud} = 1/2*(- c^{\dag}_{j\up}c_{i\up} c^{\dag}_{i\dn}c_{j\dn}
*                 - n_{i\up} n_{j_dn})
* P_{du} is obtained from P_{ud} by interchanging spin indices. 
*-----------------------------------------------------------------------*/

  // ni_nj term
  double ninj_term;
  if (op_ni_up(i)==1 && op_ni_dn(j)==1)
    ninj_term = -0.5;
  else if (op_ni_dn(i)==1 && op_ni_up(j)==1)
    ninj_term = -0.5;
  else ninj_term = 0.0;
  if (i == j) return amplitude_t(ninj_term);
  if (!op_exchange_ud(i,j)) return amplitude_t(ninj_term);

  int upspin = which_upspin();
  int up_tosite = which_upsite();
  int dnspin = which_dnspin();
  int dn_tosite = which_dnsite();
  // det_ratio for the term
  wf.get_amplitudes(psi_row, up_tosite, dnspin_sites());
  amplitude_t det_ratio1 = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  // now for dnspin hop 
  wf.get_amplitudes(psi_col, upspin_sites(), dn_tosite);
  // since the upspin should have moved
  wf.get_amplitudes(psi_col(upspin), up_tosite, dn_tosite);
  // updated 'dnspin'-th row of psi_inv
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;

  // for safety: if 'det_ratio1 == 0', result is zero
  if (std::isinf(std::abs(ratio_inv))) {
    return amplitude_t(ninj_term);
  }

  // elements other than 'upspin'-th
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
  }
  for (int i=upspin+1; i<num_upspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
  }
  inv_row(upspin) = ratio_inv * psi_inv(dnspin,upspin);
  // ratio for the dnspin hop
  amplitude_t det_ratio2 = psi_col.cwiseProduct(inv_row).sum();
  amplitude_t det_ratio = ampl_part(std::conj(det_ratio1*det_ratio2));
  /*
  if (std::isnan(det_ratio)) {
    std::cout << std::scientific<< det_ratio1 << "\n\n";
    std::cout << std::scientific<< ratio_inv << "\n\n";
    std::cout << std::scientific<< det_ratio2 << "\n\n";
    std::cout << "NaN detected\n"; getchar();
  }*/
  return -0.5 * det_ratio + amplitude_t(ninj_term);
}

amplitude_t SysConfig::apply_bondsinglet_hop(const int& fr_site_i, 
  const int& fr_site_ia, const int& to_site_j, const int& to_site_jb) const
{
/*----------------------------------------------------------------------
* Evaluates the following operator:
* Denoting: 
*          c^{dag}_{i\sigma} = d_{i\sigma}
*          c_{i\sigma}       = c_{i\sigma}
*          1/\sqrt{2}        = 1sqrt2 
*          1/2               = half 
*
* F_{ab}(i,j) = 1sqrt2 x (d_{i\up}d_{i+a\dn} - d_{i\dn}d_{i+a\up}) x
*               1sqrt2 x (c_{i\dn}c_{i+b\up} - c_{i\up}c_{i+b\dn}) 
*
* It can be written as sum of 4 terms:
* 
* F_{ab}(i,j) = 0.5 x [ (d_{i+a\dn}c_{j+b\dn) (d_{i\up}c_{j\up}) 
*                     + (d_{i+a\dn}c_{j\dn})  (d_{i\up}c_{j+b\up})
*                     + (d_{i\dn}c_{j+b\dn})  (d_{i+a\up}c_{j\up}) 
*                     + (d_{i\dn}c_{j\dn})    (d_{i+a\up}c_{j+b\up}) ]
*
* Each term describes a 'UP'-spin hop followed by a 'DOWN'-spin hop
*
* Assumption: No hopping crosses boundary (no need to consider BC twists)
*-----------------------------------------------------------------------*/
  int num_terms = 4;
  int up_frsite, up_tosite;
  int dn_frsite, dn_tosite;
  //int phase = 1;
  amplitude_t net_ratio(0.0);
  for (int iterm=0; iterm<num_terms; ++iterm) {
    switch (iterm) {
      case 0:
        up_frsite = fr_site_i;  up_tosite = to_site_j;
        dn_frsite = fr_site_ia; dn_tosite = to_site_jb;
        break;
      case 1:
        up_frsite = fr_site_i;  up_tosite = to_site_jb;
        dn_frsite = fr_site_ia; dn_tosite = to_site_j;
        break;
      case 2:
        up_frsite = fr_site_ia; up_tosite = to_site_j;
        dn_frsite = fr_site_i;  dn_tosite = to_site_jb;
        break;
      case 3:
        up_frsite = fr_site_ia; up_tosite = to_site_jb;
        dn_frsite = fr_site_i;  dn_tosite = to_site_j;
        break;
    }
/*----------------------------------------------------------------------
* Assumptions (true for bond-single hop):                     
*   up_frsite != dn_frsite
*   up_tosite != dn_tosite
*-----------------------------------------------------------------------*/
    if (up_frsite == dn_frsite) {
      throw std::range_error("SysConfig::apply_bondsinglet_hop: violation of assumption-1");
    }
    if (up_tosite == dn_tosite) {
      throw std::range_error("SysConfig::apply_bondsinglet_hop: violation of assumption-2");
    }

    int upspin, dnspin;
    amplitude_t det_ratio1, det_ratio2, term_ratio;

    // if same site for upspin
    if (up_frsite == up_tosite) {
      if (op_ni_up(up_frsite) == 0) continue;
      // if same site for dnspin
      if (dn_frsite == dn_tosite) {
        if (op_ni_dn(dn_frsite) == 0) continue;
        // term ratio
        term_ratio = amplitude_t(1.0);
      }
      else {
        // if dnspin hop to different site
        if (!op_cdagc_dn(dn_frsite,dn_tosite)) continue;
        dnspin = which_dnspin();
        // det_ratio for the term
        wf.get_amplitudes(psi_col,upspin_sites(),dn_tosite);
        det_ratio2 = psi_col.cwiseProduct(psi_inv.row(dnspin)).sum();
        // term ratio
        term_ratio = det_ratio2 * pj.gw_ratio(*this, dn_frsite, dn_tosite);
      }
    }
    else {
      // upspin hop to different site
      if (!op_cdagc_up(up_frsite, up_tosite)) continue;
      upspin = which_upspin();
      wf.get_amplitudes(psi_row, up_tosite, dnspin_sites());
      det_ratio1 = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
      // term ratio
      term_ratio = det_ratio1 * pj.gw_ratio(*this, up_frsite, up_tosite);

      // same site for dnspin
      if (dn_frsite == dn_tosite) {
        if (op_ni_dn(dn_frsite) == 0) continue;
        // above term_ratio does not change
      }
      else {
        // dnspin hop to different site
        if (std::abs(det_ratio1) < dratio_cutoff()) continue; // for safety
        if (!op_cdagc_dn(dn_frsite, dn_tosite)) continue;
        dnspin = which_dnspin();
        wf.get_amplitudes(psi_col, upspin_sites(), dn_tosite);
        // since one upspin have moved
        wf.get_amplitudes(psi_col(upspin), up_tosite, dn_tosite);
        // updated 'dnspin'-th row of psi_inv
        amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;
        // elements other than 'upspin'-th
        for (int i=0; i<upspin; ++i) {
          amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
          inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
        }
        for (int i=upspin+1; i<num_upspins_; ++i) {
          amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
          inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
        }
        inv_row(upspin) = ratio_inv * psi_inv(dnspin,upspin);
        // ratio for the dnspin hop
        det_ratio2 = psi_col.cwiseProduct(inv_row).sum();

        // term ratio (combined with the ratio for upspin hop)
        term_ratio *= det_ratio2 * pj.gw_ratio(*this, dn_frsite, dn_tosite);
      }
    }
    // net ratio due to all non-zero terms
    net_ratio += ampl_part(std::conj(term_ratio));
  }
  //int bc_phase = bphase_i * bphase_j;
  //return 0.5 * bc_phase * net_ratio;
  return 0.5 * net_ratio;
}

amplitude_t SysConfig::apply_sitepair_hop(const int& fr_site, const int& to_site) const
{
/*----------------------------------------------------------------------
* Evaluates the following operator.
* Denote: 
*          c^{dag}_{i\sigma} = d_{i\sigma}
*          c_{i\sigma}       = c_{i\sigma}
*
* F_{ab}(i,j) = (d_{i\up}d_{i\dn}) (c_{j\dn}c_{j\up}
*             = (d_{i\dn}c_{j\dn}) (d_{i\up}c_{j\up)
*
* The term describes a 'UP'-spin hop followed by a 'DOWN'-spin hop 
* between same sites
*-----------------------------------------------------------------------*/
  if (fr_site == to_site) {
    return ampl_part(op_ni_updn(fr_site));
  }

  int upspin, dnspin;
  amplitude_t det_ratio1, det_ratio2;
  amplitude_t net_ratio(0.0);

  // upspin hop
  if (!op_cdagc_up(fr_site, to_site)) return net_ratio;
  upspin = which_upspin();
  wf.get_amplitudes(psi_row, to_site, dnspin_sites());
  det_ratio1 = psi_row.cwiseProduct(psi_inv.col(upspin)).sum();
  if (std::abs(det_ratio1) < dratio_cutoff()) return net_ratio; // for safety

  // dnspin hop
  if (!op_cdagc_dn(fr_site, to_site)) return net_ratio;
  dnspin = which_dnspin();
  wf.get_amplitudes(psi_col, upspin_sites(), to_site);
  // since one upspin have moved
  wf.get_amplitudes(psi_col(upspin), to_site, to_site);
  // updated 'dnspin'-th row of psi_inv
  amplitude_t ratio_inv = amplitude_t(1.0)/det_ratio1;
  // elements other than 'upspin'-th
  for (int i=0; i<upspin; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
  }
  for (int i=upspin+1; i<num_upspins_; ++i) {
    amplitude_t beta = ratio_inv*psi_row.cwiseProduct(psi_inv.col(i)).sum();
    inv_row(i) = psi_inv(dnspin,i) - beta * psi_inv(dnspin,upspin);
  }
  inv_row(upspin) = ratio_inv * psi_inv(dnspin,upspin);
  // ratio for the dnspin hop
  det_ratio2 = psi_col.cwiseProduct(inv_row).sum();

  // net ratio for up & dn spin hop
  double gw_ratio = pj.gw_ratio_pairhop(fr_site, to_site);
  net_ratio = ampl_part(std::conj(det_ratio1*det_ratio2)*gw_ratio);

  return net_ratio;
}

void SysConfig::get_grad_logpsi(RealVector& grad_logpsi) const
{
  // grad_logpsi wrt pj parameters
  int p = pj.varparms().size();
  if (p > 0) {
    RealVector grad(p);
    pj.get_grad_logp(*this, grad);
    for (int n=0; n<p; ++n) {
      grad_logpsi(n) = grad(n);
    }
  }
  /*
  for (int n=0; n<p; ++n) {
    if (pj.varparms()[n].name()=="gfactor") {
      double g = pj.varparms()[n].value();
      grad_logpsi(n) = static_cast<double>(dblocc_count())/g;
    }
    else {
      throw std::range_error("SysConfig::get_grad_logpsi: this pj parameter not implemented\n");
    }
  }
  */

  // grad_logpsi wrt wf parameters
  for (int n=0; n<wf.varparms().size(); ++n) {
    wf.get_gradients(psi_grad,n,upspin_sites(),dnspin_sites());
    grad_logpsi(p+n) = std::real(psi_grad.cwiseProduct(psi_inv.transpose()).sum());
  }
}

double SysConfig::accept_ratio(void)
{
  // acceptance ratio wrt particle number
  return static_cast<double>(last_accepted_moves_)/
         static_cast<double>(num_upspins_+num_dnspins_); 
  //return static_cast<double>(last_accepted_moves_)/
  //       static_cast<double>(last_proposed_moves_); 
}

void SysConfig::reset_accept_ratio(void)
{
  last_proposed_moves_ = 0;
  last_accepted_moves_ = 0;
}

void SysConfig::print_stats(std::ostream& os) const
{
  long proposed_hops = num_proposed_moves_[move_t::uphop]
                         + num_proposed_moves_[move_t::dnhop];
  long proposed_exch = num_proposed_moves_[move_t::exch];
  long accepted_hops = num_accepted_moves_[move_t::uphop] 
                     + num_accepted_moves_[move_t::dnhop];
  long accepted_exch = num_accepted_moves_[move_t::exch];
  double accept_ratio, hop_ratio, exch_ratio;
  if (proposed_hops+proposed_exch>0) {
    accept_ratio = 100.0*double(accepted_hops+accepted_exch)/(proposed_hops+proposed_exch);
  }
  else  accept_ratio = 0.0;
  if (proposed_hops>0) {
    hop_ratio = double(100.0*accepted_hops)/(proposed_hops);
  }
  else hop_ratio = 0.0;
  if (proposed_exch>0) {
    exch_ratio = double(100.0*accepted_exch)/(proposed_exch);
  }
  else exch_ratio = 0.0;
  std::ios state(NULL);
  state.copyfmt(os);
  os <<std::fixed<<std::setprecision(1)<<std::right;
  os << "--------------------------------------\n";
  os << " total mcsteps = " << num_updates_ <<"\n";
  os << " total accepted moves = " << (accepted_hops+accepted_exch)<<"\n";
  os << " acceptance ratio = " << accept_ratio << " %\n";
  os << " hopping = " << hop_ratio << " %\n";
  os << " exchange = " << exch_ratio << " %\n";
  os << "--------------------------------------\n";
  os.copyfmt(state);
}


} // end namespace vmc



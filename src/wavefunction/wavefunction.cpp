/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-01-30 18:54:09
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-11 10:28:20
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "wavefunction.h"
#include <boost/algorithm/string.hpp>

namespace var {

Wavefunction::Wavefunction(const lattice::Lattice& lattice,
  const input::Parameters& inputs, const model::Hamiltonian& model, const bool& site_disorder)
  : num_sites_(lattice.num_sites())
{
  name_ = inputs.set_value("wavefunction", "NONE");
  boost::to_upper(name_);
  using order_t = MF_Order::order_t;
  using pairing_t = MF_Order::pairing_t;
  pairwf_ = true;
  if (name_ == "FERMISEA") {
    groundstate_.reset(new Fermisea(order_t::null,inputs,lattice,model));
    pairwf_ = false;
  }
  else if (name_ == "AF") {
    groundstate_.reset(new Fermisea(order_t::AF,inputs,lattice,model));
    pairwf_ = false;
  }
  else if (name_ == "SC_SWAVE") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::SWAVE,inputs,lattice,model));
  }
  else if (name_ == "SC_EXTENDED_S") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::EXTENDED_S,inputs,lattice,model));
  }
  else if (name_ == "SC_DWAVE") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::DWAVE,inputs,lattice,model));
  }
  else if (name_ == "SC_DXY") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::DXY,inputs,lattice,model));
  }
  else if (name_ == "SC_D+ID") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::D_PLUS_ID,inputs,lattice,model));
  }
  else if (name_ == "CUSTOM_SC") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::CUSTOM,inputs,lattice,model));
  }
  else if (name_ == "SC_CDW_SDW") {
    groundstate_.reset(new BCS_State(order_t::SC,pairing_t::SC_CDW_SDW,inputs,lattice,model));
  }
  else if (name_ == "DISORDERED_SC") {
    groundstate_.reset(new DisorderedSC(order_t::SC,pairing_t::DWAVE,inputs,lattice));
  }
  else {
    throw std::range_error("Wavefunction::Wavefunction: unidefined wavefunction");
  }
}

std::string Wavefunction::signature_str(void) const
{
  // signature string
  std::ostringstream signature;
  signature << "wf_N"; 
  signature << std::setfill('0'); 
  signature << std::setw(3) << groundstate_->num_upspins(); 
  signature << std::setw(3) << groundstate_->num_dnspins(); 
  return signature.str();
}

void Wavefunction::get_vparm_names(std::vector<std::string>& vparm_names, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_names[start_pos+i] = p.name(); ++i;
  }
}

void Wavefunction::get_vparm_values(var::parm_vector& vparm_values, 
  unsigned start_pos)
{
  unsigned i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void Wavefunction::get_vparm_vector(std::vector<double>& vparm_values, 
  unsigned start_pos)
{
  unsigned i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void Wavefunction::get_vparm_lbound(var::parm_vector& vparm_lb, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_lb[start_pos+i] = p.lbound(); ++i;
  }
}

void Wavefunction::get_vparm_ubound(var::parm_vector& vparm_ub, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : groundstate_->varparms()) {
    vparm_ub[start_pos+i] = p.ubound(); ++i;
  }
}


/*---------------------------------------------------------------------
*  Computation of waefunction amplitudes 
*----------------------------------------------------------------------*/
int Wavefunction::compute(const lattice::Lattice& lattice, 
  const input::Parameters& inputs, const bool& psi_gradient)
{
  groundstate_->update(inputs);
  compute_amplitudes(psi_gradient);
  return 0;
}

int Wavefunction::compute(const lattice::Lattice& lattice, const var::parm_vector& pvector,
  const unsigned& start_pos, const bool& psi_gradient)
{
  groundstate_->update(pvector,start_pos);
  compute_amplitudes(psi_gradient);
  return 0;
}

// recompute for change in lattice BC 
int Wavefunction::recompute(const lattice::Lattice& lattice)
{
  //std::cout << "recomputing\n"; 
  groundstate_->update(lattice);
  compute_amplitudes();
  /*
  std::cout << "Wavefunction::recompute\n";
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      if (std::abs(psiup_(i,j).imag())>1.0E-6) {
        std::cout << "wf["<<i<<","<<j<<"] = "<<std::imag(psiup_(i,j))<<"\n";
      }
    }
  }
  */
  return 0;
}

int Wavefunction::compute_amplitudes(const bool& psi_gradient)
{
  int num_upspins = groundstate_->num_upspins();
  int num_dnspins = groundstate_->num_dnspins();
  single_determinant_ = (num_upspins==num_dnspins);

  if (single_determinant_) {
    psiup_.resize(num_sites_,num_sites_);
    groundstate_->get_wf_amplitudes(psiup_);
    if (psi_gradient) {
      psiup_grad_.resize(varparms().size());
      for (int i=0; i<varparms().size(); ++i) {
        psiup_grad_[i].resize(num_sites_,num_sites_);
      }
      groundstate_->get_wf_gradient(psiup_grad_);
      have_gradient_ = true;
    }
    else {
      have_gradient_ = false;
    }
  }
  // product of UP and DN determinants
  else {
    psiup_.resize(num_sites_,num_upspins);
    psidn_.resize(num_dnspins,num_sites_);
    groundstate_->get_wf_amplitudes(psiup_, psidn_);
    if (psi_gradient) {
      psiup_grad_.resize(varparms().size());
      psidn_grad_.resize(varparms().size());
      for (int i=0; i<varparms().size(); ++i) {
        psiup_grad_[i].resize(num_sites_,num_upspins);
        psidn_grad_[i].resize(num_dnspins,num_sites_);
      }
      groundstate_->get_wf_gradient(psiup_grad_,psidn_grad_);
      have_gradient_ = true;
    }
    else {
      have_gradient_ = false;
    }
  }
  /*
  std::cout << "Wavefunction::compute\n";
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      if (std::abs(psiup_(i,j).imag())>1.0E-6) {
        std::cout << "psi["<<i<<","<<j<<"] = "<<std::imag(psiup_(i,j))<<"\n";
      }
    }
  }
  */
  return 0;
}

/*---------------------------------------------------------------------
*  Provides amplitude in SINGLE DETERMINANT case
*----------------------------------------------------------------------*/
void Wavefunction::get_amplitudes(Matrix& psi, const std::vector<int>& row, 
  const std::vector<int>& col) const
{
  // full matrix
  for (int i=0; i<row.size(); ++i) {
    for (int j=0; j<col.size(); ++j) {
      psi(i,j) = psiup_(row[i],col[j]);
    }
  }
}

void Wavefunction::get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const
{
  // ampl vector corresponding to a fixed 'upspin' position
  for (int j=0; j<col.size(); ++j) {
    psi_vec[j] = psiup_(irow,col[j]);
  }
}

void Wavefunction::get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const
{
  // ampl vector corresponding to a fixed 'dnspin' position
  for (int j=0; j<row.size(); ++j) {
    psi_vec[j] = psiup_(row[j],icol);
  }
}

void Wavefunction::get_amplitudes(amplitude_t& elem, const int& irow, 
  const int& jcol) const
{
  // single element corresponding to a fixed 'dnspin' & 'dnspin' position
  elem = psiup_(irow,jcol);
}

void Wavefunction::get_gradients(Matrix& psi_grad, const int& n, 
  const std::vector<int>& row, const std::vector<int>& col) const
{
  // derivative matrix
  if (!have_gradient_) {
    throw std::logic_error("Wavefunction::get_gradients: gradients were not computed");
  }
  for (int i=0; i<row.size(); ++i) {
    for (int j=0; j<col.size(); ++j) {
      psi_grad(i,j) = psiup_grad_[n](row[i],col[j]);
    }
  }
}

/*---------------------------------------------------------------------
*  Provides amplitude in PRODUCT DETERMINANT case
*----------------------------------------------------------------------*/
void Wavefunction::get_amplitudes(Matrix& psiup, Matrix& psidn, const std::vector<int>& row,  
    const std::vector<int>& col) const
{
  // full matrix corresponding all 'upspin' positions
  for (int i=0; i<row.size(); ++i) {
    psiup.row(i) = psiup_.row(row[i]);
  }
  // full matrix corresponding all 'dnspin' positions
  for (int j=0; j<col.size(); ++j) {
    psidn.col(j) = psidn_.col(col[j]);
  }
}

void Wavefunction::get_amplitudes(ColVector& psiup_vec, const int& irow) const
{
  // amplitude vector corresponding to the given 'upspin' position
  psiup_vec = psiup_.row(irow);
}

void Wavefunction::get_amplitudes(RowVector& psidn_vec, const int& icol) const
{
  // amplitude vector corresponding to the given 'dnspin' position
  psidn_vec = psidn_.col(icol);
}

void Wavefunction::get_gradients(Matrix& psiup_grad, Matrix& psidn_grad, 
  const int& n, const std::vector<int>& row, const std::vector<int>& col) const
{
  // derivative matrix
  if (!have_gradient_) {
    throw std::logic_error("Wavefunction::get_gradients: gradients were not computed");
  }
  for (int i=0; i<row.size(); ++i) {
    psiup_grad.row(i) = psiup_grad_[n].row(row[i]);
  }
  for (int j=0; j<col.size(); ++j) {
    psidn_grad.col(j) = psidn_grad_[n].col(col[j]);
  }
}


} // end namespace var












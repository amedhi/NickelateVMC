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
  if (name_ == "FERMISEA") {
    groundstate_.reset(new Fermisea(order_t::null,inputs,lattice,model));
  }
  else if (name_ == "AF") {
    groundstate_.reset(new Fermisea(order_t::AF,inputs,lattice,model));
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
  // resize
  psi_up_.resize(num_sites_,num_sites_);
  psi_gradient_.resize(varparms().size());
  for (unsigned i=0; i<varparms().size(); ++i)
    psi_gradient_[i].resize(num_sites_,num_sites_);
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

int Wavefunction::compute(const lattice::Lattice& lattice, 
  const input::Parameters& inputs, const bool& psi_gradient)
{
  groundstate_->update(inputs);
  groundstate_->get_wf_amplitudes(psi_up_);
  if (psi_gradient) {
    groundstate_->get_wf_gradient(psi_gradient_);
    have_gradient_ = true;
  }
  else have_gradient_ = false;

  /*
  std::cout << "Wavefunction::compute\n";
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      if (std::abs(psi_up_(i,j).imag())>1.0E-6) {
        std::cout << "psi["<<i<<","<<j<<"] = "<<std::imag(psi_up_(i,j))<<"\n";
      }
    }
  }
  */


  return 0;
}

int Wavefunction::compute(const lattice::Lattice& lattice, const var::parm_vector& pvector,
  const unsigned& start_pos, const bool& psi_gradient)
{
  groundstate_->update(pvector,start_pos);
  groundstate_->get_wf_amplitudes(psi_up_);
  if (psi_gradient) {
    groundstate_->get_wf_gradient(psi_gradient_);
    have_gradient_ = true;
  }
  else have_gradient_ = false;
  return 0;
}

// recompute for change in lattice BC 
int Wavefunction::recompute(const lattice::Lattice& lattice)
{
  //std::cout << "recomputing\n"; 
  groundstate_->update(lattice);
  groundstate_->get_wf_amplitudes(psi_up_);
  if (have_gradient_) {
    groundstate_->get_wf_gradient(psi_gradient_);
    //std::cout << "get_wf_gradient\n"; getchar();
  }
  /*
  std::cout << "Wavefunction::recompute\n";
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      if (std::abs(psi_up_(i,j).imag())>1.0E-6) {
        std::cout << "wf["<<i<<","<<j<<"] = "<<std::imag(psi_up_(i,j))<<"\n";
      }
    }
  }
  */
  return 0;
}

void Wavefunction::get_amplitudes(Matrix& psi, const std::vector<int>& row, 
  const std::vector<int>& col) const
{
  for (unsigned i=0; i<row.size(); ++i)
    for (unsigned j=0; j<col.size(); ++j)
      psi(i,j) = psi_up_(row[i],col[j]);
}

void Wavefunction::get_amplitudes(ColVector& psi_vec, const int& irow,  
    const std::vector<int>& col) const
{
  for (unsigned j=0; j<col.size(); ++j) {
    psi_vec[j] = psi_up_(irow,col[j]);
  }
}

void Wavefunction::get_amplitudes(RowVector& psi_vec, const std::vector<int>& row,
    const int& icol) const
{
  for (unsigned j=0; j<row.size(); ++j) {
    psi_vec[j] = psi_up_(row[j],icol);
  }
}

void Wavefunction::get_amplitudes(amplitude_t& elem, const int& irow, 
  const int& jcol) const
{
  elem = psi_up_(irow,jcol);
}

void Wavefunction::get_gradients(Matrix& psi_grad, const int& n, 
  const std::vector<int>& row, const std::vector<int>& col) const
{
  if (!have_gradient_) 
    throw std::logic_error("Wavefunction::get_gradients: gradients were not computed");
  for (unsigned i=0; i<row.size(); ++i)
    for (unsigned j=0; j<col.size(); ++j)
      psi_grad(i,j) = psi_gradient_[n](row[i],col[j]);
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

/*
int Wavefunction::compute_gradients(const lattice::Lattice& lattice, 
  const var::parm_vector& pvector, const unsigned& start_pos)
{
  // Gradient of amplitudes wrt the variational parameters 
  // by numerical differentiation (central defference formula)
  unsigned num_parm = mf_model_.varparms().size();
  psi_gradient_.resize(num_parm);
  work_mat.resize(num_sites_,num_sites_);
  //double scale = 0.005;
  unsigned i = 0;
  for (const auto& p : mf_model_.varparms()) {
    psi_gradient_[i].resize(num_sites_,num_sites_);
    //double h = scale * (p.ubound()-p.lbound());
    double h = p.diff_h();
    double inv_2h = 0.5/h;
    double x = pvector[start_pos+i];
    mf_model_.update_parameter(p.name(), x+h);
    compute_amplitudes(psi_gradient_[i], lattice);
    mf_model_.update_parameter(p.name(), x-h);
    compute_amplitudes(work_mat, lattice);
    // model to original state
    mf_model_.update_parameter(p.name(), x);
    // derivative
    psi_gradient_[i] -= work_mat;
    psi_gradient_[i] *= inv_2h;
    ++i;
  }
  return 0;
}
int Wavefunction::compute_amplitudes(Matrix& psi_mat, const lattice::Lattice& lattice)
{
  switch (type_) {
    case wf_type::bcs_oneband: 
      bcs_oneband(); 
      pair_amplitudes(lattice, psi_mat);
      break;
    case wf_type::bcs_multiband: 
      bcs_multiband(); 
      pair_amplitudes(lattice, psi_mat);
      break;
    case wf_type::bcs_disordered: 
      bcs_disordered(lattice); 
      pair_amplitudes(lattice, psi_mat);
      break;
    case wf_type::normal: 
      fermisea(); 
      fermisea_amplitudes(lattice); 
      break;
  }
  return 0;
}
void Wavefunction::pair_amplitudes(const lattice::Lattice& lattice, Matrix& psi_mat)
{
  double one_by_nk = 1.0/static_cast<double>(num_kpoints_);
  for (unsigned i=0; i<num_sites_; ++i) {
    //unsigned m = lattice.site_uid(i);
    unsigned m = blochbasis_.representative_state_idx(i);
    auto Ri = lattice.site_cellcord(i);
    for (unsigned j=0; j<num_sites_; ++j) {
      //unsigned n = lattice.site_uid(j);
      unsigned n = blochbasis_.representative_state_idx(j);
      auto Rj = lattice.site_cellcord(j);
      std::complex<double> ksum(0.0);
      for (unsigned k=0; k<num_kpoints_; ++k) {
        Vector3d kvec = blochbasis_.kvector(k);
        ksum += cphi_k[k](m,n) * std::exp(ii()*kvec.dot(Ri-Rj));
      }
      psi_mat(i,j) = ampl_part(ksum) * one_by_nk;
      //std::cout << psi_up_(i,j) << "\n"; 
      //getchar();
    }
  }
}
*/



} // end namespace var












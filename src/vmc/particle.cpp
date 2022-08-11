/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2022-08-11 12:49:54
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./particle.h"

namespace vmc {

//-------------------Site Occupancy--------------------------------
void SiteOccupancy::setup(const lattice::Lattice& lattice, const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  num_basis_sites_ = lattice.num_basis_sites();
  num_particles_ = config.num_particles();
  std::vector<std::string> elem_names(num_basis_sites_);
  for (int i=0; i<num_basis_sites_; ++i) {
  	std::ostringstream ss;
  	ss << "site-"<<i;
  	elem_names[i] = ss.str();
  }
  this->resize(elem_names.size(), elem_names);
  //this->set_have_total();
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}

void SiteOccupancy::measure(const lattice::Lattice& lattice, const SysConfig& config) 
{
  IntVector matrix_elem(num_basis_sites_);
  IntVector num_subsites(num_basis_sites_);
  matrix_elem.setZero();
  num_subsites.setZero();
  //for (auto s=graph.sites_begin(); s!=graph.sites_end(); ++s) {
  //  int site = graph.site(s);
  //  int basis = graph.site_uid(s);
  for (const auto& s : lattice.sites()) {
    int site = s.id();
    int basis = s.uid();
    matrix_elem(basis) += config.apply(model::op::ni_sigma(),site);
    num_subsites(basis) += 1;
  }
  for (int i=0; i<num_basis_sites_; ++i) {
  	//config_value_[i] = static_cast<double>(matrix_elem[i])/num_particles_;
    config_value_[i] = static_cast<double>(matrix_elem[i])/num_subsites[i];
  }
  // add to databin
  *this << config_value_;
}




} // end namespave vmc









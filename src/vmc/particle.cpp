/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2022-11-12 23:10:52
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./particle.h"

namespace vmc {

//-------------------Site Particle Density--------------------------------
void ParticleDensity::setup(const lattice::Lattice& lattice, const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  num_site_types_ = lattice.num_site_types();
  num_particles_ = config.num_particles();
  std::vector<std::string> elem_names(num_site_types_+1);
  elem_names[0] = "total";
  for (int i=0; i<num_site_types_; ++i) {
    std::ostringstream ss;
    ss << "<n>["<<i<<"]";
    elem_names[i+1] = ss.str();
  }
  this->resize(elem_names.size(), elem_names);
  //this->set_have_total();
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}

void ParticleDensity::measure(const lattice::Lattice& lattice, const SysConfig& config) 
{
  IntVector matrix_elem(num_site_types_);
  IntVector num_subsites(num_site_types_);
  matrix_elem.setZero();
  num_subsites.setZero();
  for (const auto& s : lattice.sites()) {
    int site = s.id();
    int type = s.type();
    matrix_elem(type) += config.apply(model::op::ni_sigma(),site);
    num_subsites(type) += 1;
  }
  double total = 0.0;
  for (int i=0; i<num_site_types_; ++i) {
    config_value_[i+1] = static_cast<double>(matrix_elem[i])/num_subsites[i];
    total += config_value_[i+1];
  }
  config_value_[0] = total;

  // add to databin
  *this << config_value_;
}


//-------------------Doublon & Holon Density--------------------------------
void DoublonDensity::setup(const lattice::Lattice& lattice, const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  num_site_types_ = lattice.num_site_types();
  std::vector<std::string> elem_names(num_site_types_*2);
  for (int i=0; i<num_site_types_; ++i) {
    std::ostringstream ss;
    ss << "<dblon>["<<i<<"]";
    elem_names[2*i] = ss.str();
    ss.clear(); ss.str("");
    ss << "<holon>["<<i<<"]";
    elem_names[2*i+1] = ss.str();
  }
  this->resize(elem_names.size(), elem_names);
  //this->set_have_total();
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}

void DoublonDensity::measure(const lattice::Lattice& lattice, const SysConfig& config) 
{
  IntVector dblon_count(num_site_types_);
  IntVector holon_count(num_site_types_);
  IntVector num_subsites(num_site_types_);
  dblon_count.setZero();
  holon_count.setZero();
  num_subsites.setZero();
  for (const auto& s : lattice.sites()) {
    int site = s.id();
    int type = s.type();
    dblon_count[type] += config.op_ni_dblon(site);
    holon_count[type] += config.op_ni_holon(site);
    num_subsites[type] += 1;
  }
  for (int i=0; i<num_site_types_; ++i) {
    config_value_[2*i] = static_cast<double>(dblon_count[i])/num_subsites[i];
    config_value_[2*i+1] = static_cast<double>(holon_count[i])/num_subsites[i];
  }

  // add to databin
  *this << config_value_;
}


} // end namespave vmc









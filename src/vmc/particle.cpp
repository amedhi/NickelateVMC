/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2023-04-22 23:16:41
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
  std::vector<std::string> elem_names; 
  elem_names.push_back("total");
  for (int i=0; i<num_site_types_; ++i) {
    std::ostringstream ss;
    ss << "<n>["<<i<<"]";
    elem_names.push_back(ss.str());
  }
  // Sz value
  for (int i=0; i<num_site_types_; ++i) {
    std::ostringstream ss;
    ss << "<Sz>["<<i<<"]";
    elem_names.push_back(ss.str());
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
  IntVector Sz(num_site_types_);
  Sz.setZero();
  for (const auto& s : lattice.sites()) {
    int site = s.id();
    int type = s.type();
    Sz(type) += config.apply(model::op::Sz(),site);
  }

  double total = 0.0;
  for (int i=0; i<num_site_types_; ++i) {
    config_value_[i+1] = static_cast<double>(matrix_elem[i])/num_subsites[i];
    total += config_value_[i+1];
  }
  config_value_[0] = total;

  // Sz avg
  for (int i=0; i<num_site_types_; ++i) {
    config_value_[num_site_types_+1+i] = static_cast<double>(Sz[i])/num_subsites[i];
  }

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
    dblon_count[type] += config.apply_ni_dblon(site);
    holon_count[type] += config.apply_ni_holon(site);
    num_subsites[type] += 1;
  }
  for (int i=0; i<num_site_types_; ++i) {
    config_value_[2*i] = static_cast<double>(dblon_count[i])/num_subsites[i];
    config_value_[2*i+1] = static_cast<double>(holon_count[i])/num_subsites[i];
  }

  // add to databin
  *this << config_value_;
}

//-------------------Momentum Distribution--------------------------------
void MomentumDist::setup(const lattice::Lattice& lattice, const SysConfig& config)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  num_site_types_ = lattice.num_site_types();
  num_particles_ = config.num_particles();
  std::vector<std::string> elem_names; 
  elem_names.push_back("total");
  for (int i=0; i<num_site_types_; ++i) {
    std::ostringstream ss;
    ss << "<nk>["<<i<<"]";
    elem_names.push_back(ss.str());
  }
  this->resize(elem_names.size(), elem_names);
  //this->set_have_total();
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}

void MomentumDist::measure(const lattice::Lattice& lattice, const SysConfig& config) 
{
  config_value_.setZero();

  // add to databin
  *this << config_value_;
}


} // end namespave vmc









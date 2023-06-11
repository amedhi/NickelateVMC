/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2023-06-11 23:32:21
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
  // store kpoints
  num_kpoints_ = config.wavefunc().blochbasis().num_kpoints();
  kpoints_.resize(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    kpoints_[k] = config.wavefunc().blochbasis().kvector(k);
  }
  nk_data_.resize(num_kpoints_,num_site_types_);

  std::vector<std::string> elem_names; 
  //elem_names.push_back("total");
  for (int i=0; i<num_site_types_; ++i) {
    std::ostringstream ss;
    ss << "<nk>["<<i<<"]";
    elem_names.push_back(ss.str());
  }
  this->resize(nk_data_.size(), elem_names);
  //config_value_.resize(elem_names.size());
  setup_done_ = true;
}

void MomentumDist::measure(const lattice::Lattice& lattice, const SysConfig& config) 
{
  //config_value_.setZero();
  nk_data_.setZero();
  std::complex<double> ampl;

  // diagonal terms
  for (const auto& s : lattice.sites()) {
    int ni = config.apply(model::op::ni_sigma(),s.id());
    for (unsigned k=0; k<num_kpoints_; ++k) {
      nk_data_(k,s.type()) += 0.5*ni;
    }
  }

  // off-diagonal terms
  for (const auto& s1 : lattice.sites()) {
    int stype = s1.type();
    Vector3d R1 = s1.coord();
    for (const auto& s2 : lattice.sites()) {
      if (s1.id() == s2.id()) continue;
      if (s2.type() != stype) continue;
      Vector3d R = R1 - s2.coord();
      //std::cout << s1.id() << " - " << s2.id() << "\n";
      //std::cout << s1.type() << " - " << s2.type() << "\n\n";
      //getchar();

      // ampl = 0.5*\sum_{\sigma} <c^\dag_{ia\sigma}c_{ja\sigma}> 
      ampl = config.apply_cdagc_up(s1.id(),s2.id(),1,1.0);
      ampl += config.apply_cdagc_dn(s1.id(),s2.id(),1,1.0);
      ampl *= 0.5;
      for (unsigned k=0; k<num_kpoints_; ++k) {
        Vector3d kvec = kpoints_[k];
        nk_data_(k,stype) += std::real(std::exp(ii()*kvec.dot(R)) * ampl);
      }
    }
  }

  // normalization
  nk_data_ /= num_kpoints_;

  // reshape add to databin
  *this << Eigen::Map<mcdata::data_t>(nk_data_.data(), nk_data_.size());
}

void MomentumDist::print_heading(const std::string& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!is_open()) open_file();
  fs_ << header;
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;

  fs_ << std::setw(15)<<"kx";
  fs_ << std::setw(15)<<"ky";
  fs_ << std::setw(15)<<"kz";

  for (const auto& name : elem_names_) fs_ << std::setw(15)<<name<<std::setw(11)<<"err";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "#  ";
  for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  fs_ << std::endl;
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void MomentumDist::print_result(const std::vector<double>& xvals) 
{
  if (!is_on()) return;
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  fs_ << "#";
  for (const auto& p : xvals) fs_ << std::setw(14) << p;
  fs_ << "\n#" << std::string(72, '-') << "\n";

  for (int k=0; k<num_kpoints_; ++k) {
    //for (const auto& p : xvals) fs_ << std::setw(14) << p;
    fs_ << std::scientific << std::setw(15) << kpoints_[k](0);
    fs_ << std::scientific << std::setw(15) << kpoints_[k](1);
    fs_ << std::scientific << std::setw(15) << kpoints_[k](2);
    int n = k;
    for (int p=0; p<num_site_types_; ++p) {
      //fs_ << MC_Data::result_str(n);
      fs_ << std::scientific << std::setw(15) << std::abs(MC_Data::mean(n));
      fs_ << std::fixed << std::setw(10) << MC_Data::stddev(n);
      n += num_kpoints_;
    }
    fs_ << MC_Data::conv_str(k); 
    fs_ << std::endl; 
    if (k+1<num_kpoints_) {
      if (std::abs(kpoints_[k](1)-kpoints_[k+1](1)) > 1.0E-6) {
        fs_ << std::endl;
      }
    }
  }

  fs_ << std::flush;
  close_file();
}


} // end namespave vmc









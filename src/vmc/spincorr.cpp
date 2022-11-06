/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:35:10
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:08:00
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./spincorr.h"
#include "../utils/curve_fit.h"

namespace vmc {

void SpinCorrelation::setup(const lattice::Lattice& lattice)
{
  // max distance (catch sites along )
  Vector3d R = lattice.size1() * lattice.vector_a1();
  max_dist_ = std::nearbyint(R.norm());
  sitepair_list_.resize(max_dist_+1);

  for (int i=0; i<lattice.num_sites(); ++i) {
    Vector3i n1 = lattice.site(i).bravindex();
    Vector3d r1 = lattice.site(i).coord();
    for (int j=i; j<lattice.num_sites(); ++j) {
      Vector3i n2 = lattice.site(j).bravindex();
      if (n1(1)!=n2(1) || n1(2)!=n2(2)) break;
      Vector3d R = r1-lattice.site(j).coord();
      if (std::abs(R(1))>1.0E-6 or std::abs(R(2))>1.0E-6) continue;
      int d = std::nearbyint(std::abs(R(0)));
      if (d > max_dist_) continue;
      sitepair_list_.pairs_at_dist(d).push_back({i,j});
      //std::cout<<i<<" -- "<<j<<" = "<<d<<"\n"; //getchar();
    }
  }

  //corr_data_.resize(bondpair_types_.size(), max_dist_);
  corr_data_.resize(max_dist_+1);
  count_.resize(max_dist_+1);
  std::vector<std::string> elem_names(1);
  elem_names[0] = "<Sz_i.Sz_j>";
  this->resize(corr_data_.size(), elem_names);
}

void SpinCorrelation::measure(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  corr_data_.setZero();
  count_.setZero();
  for (int d=0; d<=max_dist_; ++d) {
    for (const auto& p : sitepair_list_.pairs_at_dist(d)) {
      int site_i = p.first;
      int site_j = p.second;
      int Sz_i = config.apply(model::op::Sz(),site_i);
      int Sz_j = config.apply(model::op::Sz(),site_j);
      corr_data_(d) += Sz_i*Sz_j;
      count_(d) += 1;
    }
  }

  for (int d=0; d<=max_dist_; ++d) {
    corr_data_(d) *= 0.5;
    if (count_(d) != 0) corr_data_(d) /= count_(d);
  }
  corr_data_(max_dist_) = corr_data_(0);

  // reshape add to databin
  *this << corr_data_;
}


void SpinCorrelation::print_heading(const std::string& header,
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
  // total value
  for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  fs_ << std::setw(6)<<"d";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  /*if (xvars.size()>0) {
    fs_ << "# ";
    for (const auto& p : xvars) fs_<<std::setw(14)<<p.substr(0,14);
    fs_ << std::endl;
  }*/
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void SpinCorrelation::print_result(const std::vector<double>& xvals) 
{
  if (!is_on()) return;
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);

  for (int d=0; d<=max_dist_; ++d) {
    for (const auto& p : xvals) fs_ << std::setw(14) << p;
    fs_ << std::setw(6) << d; 
    fs_ << MC_Data::result_str(d);
    fs_ << MC_Data::conv_str(d); 
    fs_ << std::endl; 
  }

  // structure factor
  double S = 0.0;
  int sign = 1;
  for (int d=0; d<=max_dist_; ++d) {
    S += sign*MC_Data::mean_data()(d);
    sign = -sign;
  }

  fs_ << "# Structure factor S(pi,pi)\n";
  fs_ << "# ";
  for (const auto& p : xvals) fs_ << std::setw(14) << p;
  fs_ << std::setw(14)<<S<<"\n";

  //fstream() << MC_Data::conv_str(0); //.substr(0,10); 
  fs_ << std::flush;
  close_file();
}


} // end namespace vmc


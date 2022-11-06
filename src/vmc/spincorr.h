/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:32:31
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:07:57
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_SPINCORR_H
#define OBS_SPINCORR_H

#include <string>
#include <vector>
#include <stdexcept>
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "../mcdata/mc_observable.h"
#include "./sysconfig.h"

namespace vmc {

class SitePairs
{
public:
  using site_t = int;
  using site_pair_t = std::pair<site_t,site_t>;
  using site_pair_list = std::vector<site_pair_t>;
  SitePairs() { data_.clear(); }
  ~SitePairs() {}
  void resize(const int& max_dist) { data_.resize(max_dist); }
  void clear(void) { data_.clear(); }
  site_pair_list& pairs_at_dist(const int& dist) { return data_[dist]; }
  const site_pair_list& pairs_at_dist(const int& dist) const { return data_[dist]; }
private:
  std::vector<site_pair_list> data_; 
};

class SpinCorrelation : public mcdata::MC_Observable
{
public:
  using site_t = unsigned;
  using MC_Observable::MC_Observable;
  void setup(const lattice::Lattice& lattice);
  void measure(const lattice::Lattice& lattice, const model::Hamiltonian& model,
    const SysConfig& config);
  void print_heading(const std::string& header, 
    const std::vector<std::string>& xvars) override;
  void print_result(const std::vector<double>& xvals) override; 
private:
  int max_dist_{0};
  SitePairs sitepair_list_;
  Eigen::VectorXd corr_data_;
  Eigen::VectorXi count_;
};



} // end namespace vmc

#endif

/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:32:31
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:07:57
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_SCCORR_H
#define OBS_SCCORR_H

#include <string>
#include <vector>
#include <stdexcept>
#include "../lattice/graph.h"
#include "../model/model.h"
#include "../mcdata/mc_observable.h"
#include "./sysconfig.h"

namespace vmc {

class symmetry_site_pairs
{
public:
  using site_t = lattice::LatticeGraph::site_descriptor;
  using site_pair_t = std::pair<site_t,site_t>;
  using site_pair_list = std::vector<site_pair_t>;
  symmetry_site_pairs() { data_.clear(); }
  ~symmetry_site_pairs() {}
  void resize(const int& max_dist) { data_.resize(max_dist); }
  void clear(void) { data_.clear(); }
  site_pair_list& pairs_at_dist(const int& dist) { return data_[dist]; }
  const site_pair_list& pairs_at_dist(const int& dist) const { return data_[dist]; }
private:
  std::vector<site_pair_list> data_; 
};

class SC_Correlation : public mcdata::MC_Observable
{
public:
  using site_t = lattice::LatticeGraph::site_descriptor;
  using MC_Observable::MC_Observable;
  void setup(const lattice::LatticeGraph& graph, const var::MF_Order::pairing_t& pairsymm);
  void measure(const lattice::LatticeGraph& graph, const model::Hamiltonian& model,
    const SysConfig& config);
  //const unsigned& num_site_pairs(void) const { return src_pairs_size_; }
  //const std::pair<site_t,site_t>& site_pair(const unsigned& i) const 
  //  { return src_pairs_[i]; }
  void print_heading(const std::string& header, 
    const std::vector<std::string>& xvars) override;
  void print_result(const std::vector<double>& xvals) override; 
private:
  var::MF_Order::pairing_t pair_symm_{var::MF_Order::pairing_t::DWAVE};
  int num_basis_sites_{1};
  int max_dist_{0};
  symmetry_site_pairs symm1_pairs_;
  symmetry_site_pairs symm2_pairs_;
  std::vector<symmetry_site_pairs> symm_list_;

  std::vector<std::pair<int,int> > bondpair_types_;
  Eigen::MatrixXd corr_data_;

  /*
  unsigned num_bond_types_{0};
  unsigned src_pairs_size_{0};
  std::vector<std::pair<site_t,site_t> > src_pairs_;
  std::vector<int> pair_distance_;
  std::vector<Eigen::MatrixXd> bond_pair_corr_;
  std::vector<Eigen::MatrixXi> num_symm_pairs_;
  */
  void measure_swave(const lattice::LatticeGraph& graph, const model::Hamiltonian& model,
    const SysConfig& config);
  void measure_dwave(const lattice::LatticeGraph& graph, const model::Hamiltonian& model,
    const SysConfig& config);
};

// to fit F(r)
class FitFunc
{
public:
  FitFunc(void) {}
  ~FitFunc(void) {}
  double operator()(const double& x, const Eigen::VectorXd& p)
  {
    return p[0] + p[1]*std::exp(-x/p[2])*std::cos(p[3]*x);
  }
  void derivative(const double& x, const Eigen::VectorXd& p, Eigen::VectorXd& pgrad)
  {
    double xp = std::exp(-x/p[2]);
    double cs = std::cos(p[3]*x);
    pgrad[0] = 1.0;
    pgrad[1] = xp*cs; 
    pgrad[2] = p[1]*pgrad[1]*x/(p[2]*p[2]);
    pgrad[3] = -p[1]*xp*std::sin(p[3]*x)*x;
  }
};

class FitFunc2
{
public:
  FitFunc2(void) {}
  ~FitFunc2(void) {}
  double operator()(const double& x, const Eigen::VectorXd& p)
  {
    return p[0] + p[1]*std::exp(-x/p[2])*std::cos(x);
  }
  void derivative(const double& x, const Eigen::VectorXd& p, Eigen::VectorXd& pgrad)
  {
    double xp = std::exp(-x/p[2]);
    double cs = std::cos(x);
    pgrad[0] = 1.0;
    pgrad[1] = xp*cs; 
    pgrad[2] = p[1]*pgrad[1]*x/(p[2]*p[2]);
  }
};



} // end namespace vmc

#endif

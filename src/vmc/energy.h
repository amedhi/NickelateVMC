/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:41:40
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-17 10:45:10
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_ENERGY_H
#define OBS_ENERGY_H

#include "../mcdata/mc_observable.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "./sysconfig.h"
#include "./disorder.h"

namespace vmc {

class Energy : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::Lattice& lattice, const model::Hamiltonian& model);
  void measure(const lattice::Lattice& lattice, const model::Hamiltonian& model,
    const SysConfig& config, const SiteDisorder& site_disorder);
  const mcdata::data_t& config_value(void) const { return config_value_; }
private:
  bool setup_done_{false};
  unsigned num_sites_{0};
  mcdata::data_t config_value_;
};

class EnergyGradient : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const SysConfig& config, const int& sample_size=500);
  void measure(const SysConfig& config, const double& config_energy);
  void finalize(void);
  void reset(void) override;
  void reset_batch_limit(const int& sample_size);
  const RealVector& grad_logpsi(void) const { return grad_logpsi_; }
private:
  bool setup_done_{false};
  unsigned num_varp_{0};
  int batch_size_{500};
  mcdata::data_t config_value_;
  RealVector grad_logpsi_;
  mcdata::MC_Observable batch_gsum_{"batch_gsum"};
  mcdata::MC_Observable batch_esum_{"batch_esum"};
};

class SR_Matrix : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::Lattice& lattice, const SysConfig& config,
    const int& sample_size=500);
  void reset(void) override;
  void measure(const RealVector& grad_logpsi);
  void get_matrix(Eigen::MatrixXd& sr_mat) const;
  void reset_batch_limit(const int& sample_size);
private:
  bool setup_done_{false};
  unsigned num_sites_{0};
  unsigned num_varp_{0};
  int batch_size_{500};
  mutable mcdata::data_t config_value_;
  mutable mcdata::data_t upper_elems_;
  mcdata::MC_Observable batch_sum_{"batch_sum"};
};


} // end namespace vmc

#endif

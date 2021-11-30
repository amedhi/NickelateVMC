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
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./sysconfig.h"
#include "./disorder.h"

namespace vmc {

class Energy : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::LatticeGraph& graph, const model::Hamiltonian& model);
  void measure(const lattice::LatticeGraph& graph, const model::Hamiltonian& model,
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
  void setup(const SysConfig& config);
  void measure(const SysConfig& config, const double& config_energy);
  void finalize(void);
  void reset(void) override 
  { 
    MC_Observable::reset(); grad_terms_.reset(); total_en_.reset();
  }
  const RealVector& grad_logpsi(void) const { return grad_logpsi_; }
  void MPI_send_data(const mpi::mpi_communicator& mpi_comm, const mpi::proc& proc, const int& msg_tag) override;
  void MPI_add_data(const mpi::mpi_communicator& mpi_comm, const mpi::proc& proc, const int& msg_tag) override;
private:
  bool setup_done_{false};
  unsigned num_varp_{0};
  mcdata::data_t config_value_;
  RealVector grad_logpsi_;
  mcdata::MC_Observable grad_terms_{"gradient_terms"};
  mcdata::MC_Observable total_en_{"total_energy"};
};

class SR_Matrix : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::LatticeGraph& graph, const SysConfig& config);
  void measure(const RealVector& grad_logpsi);
  void get_matrix(Eigen::MatrixXd& sr_mat) const;
private:
  bool setup_done_{false};
  unsigned num_sites_{0};
  unsigned num_varp_{0};
  mutable mcdata::data_t config_value_;
};


} // end namespace vmc

#endif

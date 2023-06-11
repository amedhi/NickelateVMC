/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-09-26 13:53:41
* @Last Modified by:   Amal Medhi, amedhi@mbpro
* @Last Modified time: 2019-09-26 13:53:41
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OBS_PARTICLE_H
#define OBS_PARTICLE_H

#include "../mcdata/mc_observable.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "./sysconfig.h"

namespace vmc {

class ParticleDensity : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::Lattice& lattice, const SysConfig& config);
  void measure(const lattice::Lattice& lattice, const SysConfig& config);
  const mcdata::data_t& config_value(void) const { return config_value_; }
private:
  bool setup_done_{false};
  int num_sites_{0};
  int num_site_types_{0};
  int num_particles_{0};
  mcdata::data_t config_value_;
};

class DoublonDensity : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::Lattice& lattice, const SysConfig& config);
  void measure(const lattice::Lattice& lattice, const SysConfig& config);
  const mcdata::data_t& config_value(void) const { return config_value_; }
private:
  bool setup_done_{false};
  int num_sites_{0};
  int num_site_types_{0};
  mcdata::data_t config_value_;
};

class MomentumDist : public mcdata::MC_Observable
{
public:
  using MC_Observable::MC_Observable;
  void setup(const lattice::Lattice& lattice, const SysConfig& config);
  void measure(const lattice::Lattice& lattice, const SysConfig& config);
  //const mcdata::data_t& config_value(void) const { return config_value_; }
  void print_heading(const std::string& header, const std::vector<std::string>& xvars) override;
  void print_result(const std::vector<double>& xvals) override; 
private:
  bool setup_done_{false};
  int num_sites_{0};
  int num_site_types_{0};
  int num_kpoints_{0};
  std::vector<Vector3d> kpoints_;
  Eigen::MatrixXd nk_data_;
  //mcdata::data_t config_value_;
};

} // end namespave vmc



#endif
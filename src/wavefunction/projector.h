/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:03:44
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-02 23:12:57
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef PROJECTOR_H
#define PROJECTOR_H

#include <iostream>
#include <vector>
#include <map>
#include "../scheduler/worker.h"
#include "../vmc/basisstate.h"
#include "../lattice/lattice.h"
#include "./varparm.h"

namespace var {

//constexpr double gw_cutoff(void) { return 1.0E-8; } 
constexpr double gw_cutoff(void) { return 1.0E-3; } 

enum class pp {gutzwiller, end};

enum class gw_ptype {DOUBLON, HOLON};

class NNTable
{
public:
  NNTable() {}
  ~NNTable() {}
  void clear(void) { nn1_list_.clear(); nn2_list_.clear(); nn3_list_.clear(); }
  void add_nn1site(const int& id) { nn1_list_.push_back(id); }
  void add_nn2site(const int& id) { nn2_list_.push_back(id); }
  void add_nn3site(const int& id) { nn3_list_.push_back(id); }
  const std::vector<int>& nn1sites(void) const { return nn1_list_; }
  const std::vector<int>& nn2sites(void) const { return nn2_list_; }
  const std::vector<int>& nn3sites(void) const { return nn3_list_; }
private:
  std::vector<int> nn1_list_;
  std::vector<int> nn2_list_;
  std::vector<int> nn3_list_;
};

class GW_Projector 
{
public:
  GW_Projector() {}
  ~GW_Projector() {}
  int init(const lattice::Lattice& lattice, const input::Parameters& inputs,
    VariationalParms& vparms);
  bool is_present(void) const { return is_present_; }
  bool is_strong(void) const; 
  bool have_holon_projection(void) const;
  int update_parameters(const VariationalParms& vparms);
  double gw_ratio(const vmc::BasisState& state, 
    const int& fr_site, const int& to_site) const;
  void get_grad_logp(const vmc::BasisState& state, RealVector& grad) const;
private:
  enum class pjn {DOUBLON, HOLON};
  bool is_present_{false};
  bool uniform_projection_{true};
  int num_sites_;
  int num_site_types_;
  std::vector<pjn> projection_; 
  std::vector<double> gw_factor_;
  std::vector<int> sitetype_list_;
  RealMatrix ratio_table_;

  void set_ratio_table(void); 
};


class WavefunProjector 
{
public:
  WavefunProjector() {}
  WavefunProjector(const lattice::Lattice& lattice, const input::Parameters& parms) 
    { init(lattice, parms); }
  ~WavefunProjector() {}
  void init(const lattice::Lattice& lattice, const input::Parameters& inputs); 
  void update(const input::Parameters& inputs); 
  void update(const var::parm_vector& pvector, const unsigned& start_pos=0);
  bool gw_projection(void) const { return gw_projector_.is_present(); }
  bool gw_projection_strong(void) const { return gw_projector_.is_strong(); }
  bool have_holon_projection(void) const { return gw_projector_.have_holon_projection(); }
  double gw_ratio(const vmc::BasisState& state, const int& fr_site, const int& to_site) const
    { return gw_projector_.gw_ratio(state, fr_site, to_site); }
  const bool& have_dh_projector(void) const { return dh_projector_; }
  double dh_factor1(void) const; 
  double dh_factor2(void) const; 
  double dh_factor3(void) const; 
  double dh_ratio(const vmc::BasisState& state, const int& fr_site, const int& to_site) const;
  void get_grad_logp(const vmc::BasisState& state, RealVector& grad) const;
  const VariationalParms& varparms(void) const { return varparms_; }
  void get_vparm_names(std::vector<std::string>& names, unsigned start_pos=0) const; 
  void get_vparm_values(var::parm_vector& values, unsigned start_pos=0) const; 
  void get_vparm_vector(std::vector<double>& vparm_values, unsigned start_pos=0) const;
  void get_vparm_lbound(var::parm_vector& lbounds, unsigned start_pos=0) const; 
  void get_vparm_ubound(var::parm_vector& ubounds, unsigned start_pos=0) const; 
private:
  using vparm_t = std::pair<std::string,double>;
  GW_Projector gw_projector_;

  bool dh_projector_{false};
  int num_sites_{0};
  int num_site_types_{1};
  int dh_range_{0};

  std::vector<double> dh_factor_;
  double dh_factor1_;
  double dh_factor2_;
  double dh_factor3_;

  //int num_gw_factors_{0};
  VariationalParms varparms_;
  std::vector<NNTable> adjacency_;

  int init_dh_projector(const lattice::Lattice& lattice, const input::Parameters& inputs); 
};


} // end namespace var

#endif
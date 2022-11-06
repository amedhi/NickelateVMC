/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:19:36
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:15:18
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef VMC_H
#define VMC_H

#include "../scheduler/worker.h"
#include "../scheduler/mpi_comm.h"
#include "../lattice/lattice.h"
#include "../model/model.h"
#include "./observables.h"
#include "./sysconfig.h"
#include "./disorder.h"
//#include "../optimizer/optimizer.h"

namespace vmc {

enum class run_mode {normal, en_function, sr_function};

class VMC //: public optimizer::Problem
{
public:
  VMC(const input::Parameters& inputs); 
  virtual ~VMC() {}
  int start(const input::Parameters& inputs, const run_mode& mode=run_mode::normal, 
    const bool& silent=false);
  int start(const var::parm_vector& varp, const run_mode& mode=run_mode::normal,
    const bool& silent=false);
  void seed_rng(const int& seed=1);
  int set_run_mode(const run_mode& mode);
  int run_simulation(void);
  int run_simulation(const int& sample_size);
  int run_simulation(const std::vector<int>& bc_list);
  int run_simulation(const int& sample_size, const std::vector<int>& bc_list);
  //int run_simulation(const Eigen::VectorXd& varp);
  int reset_observables(void);
  int do_warmup_run(void); 
  int do_measure_run(const int& num_samples); 
  int build_config(const var::parm_vector& vparms, const bool& with_psi_grad); 
  int sr_function(const var::parm_vector& vparms, double& en_mean, double& en_stddev,
    RealVector& grad, RealVector& grad_err, RealMatrix& sr_matrix,const bool& silent=false);
  int en_function(const var::parm_vector& vparms, double& en_mean, double& en_stddev,
    RealVector& grad, RealVector& grad_stddev,const bool& silent=false);
  int operator()(const var::parm_vector& vparms, double& en_mean, double& en_stddev,
     RealVector& grad, RealVector& grad_stddev) 
  { return en_function(vparms,en_mean,en_stddev,grad,grad_stddev); }
  //void get_vparm_values(var::parm_vector& varparms) 
  //  { varparms = config.vparm_values(); }
  const int& num_boundary_twists(void) const { return lattice.num_boundary_twists(); } 
  const int& num_measure_steps(void) const { return num_measure_steps_; } 
  const unsigned& num_varp(void) const { return config.num_varparms(); } 
  const var::parm_vector& varp_values(void) { return config.vparm_values(); }
  const var::parm_vector& varp_lbound(void) const { return config.vparm_lbound(); }
  const var::parm_vector& varp_ubound(void) const { return config.vparm_ubound(); }
  const std::vector<std::string>& varp_names(void) const { return config.varp_names(); }
  RandomGenerator& rng(void) const { return config.rng(); }
  const double& hole_doping(void) const { return config.hole_doping(); }
  const std::vector<std::string>& xvar_names(void) const { return xvar_names_; }
  const std::vector<double>& xvar_values(void) const { return xvar_values_; }
  const ObservableSet& observable(void) const { return observables; }
  void finalize_results(void) { observables.finalize(); } 
  void print_results(void); 
  const std::string prefix_dir(void) const { return prefix_; }
  std::ostream& print_info(std::ostream& os) const { return os << info_str_.str(); }
  static void copyright_msg(std::ostream& os);
  const bool& disordered_system(void) const { return site_disorder_.exists(); }
  const unsigned& num_disorder_configs(void) const { return site_disorder_.num_configs(); }

  void MPI_send_results(const mpi::mpi_communicator& mpi_comm, const mpi::proc& proc, 
    const int& msg_tag) { observables.MPI_send_results(mpi_comm, proc, msg_tag); }
  void MPI_recv_results(const mpi::mpi_communicator& mpi_comm, const mpi::proc& proc, 
    const int& msg_tag) { observables.MPI_recv_results(mpi_comm, proc, msg_tag); }

  // disordered case
  int disorder_start(const input::Parameters& inputs, const unsigned& disorder_config, 
    const run_mode& mode=run_mode::normal, const bool& silent=false);
  void save_optimal_parms(const var::parm_vector& optimal_parms) 
    { site_disorder_.save_optimal_parms(optimal_parms); }
  bool optimal_parms_exists(const unsigned& config) 
    { return site_disorder_.optimal_parms_exists(config); } 
  void set_disorder_config(const unsigned& config) 
    { site_disorder_.set_current_config(config); }

  // optimizer
  //void set_box_constraints(void) 
  //  { Problem::setBoxConstraint(varp_lbound(), varp_ubound()); }
private:
  run_mode run_mode_{run_mode::normal};
  lattice::Lattice lattice;
  model::Hamiltonian model;
  SysConfig config;
  SiteDisorder site_disorder_;
  int rng_seed_;
  unsigned num_sites_;
  unsigned num_varparms_;
  // default BC_TWIST list
  std::vector<int> bc_default_{0};

  // observables
  ObservableSet observables;
  std::string prefix_{"./"};
  std::vector<std::string> xvar_names_;
  std::vector<double> xvar_values_;


  // mc parameters
  enum move_t {uphop, dnhop, exch, end};
  int num_measure_steps_{0}; 
  int num_warmup_steps_{0};
  int min_interval_{0};
  int max_interval_{0};
  int check_interval_{0};
  bool silent_mode_{false};

  mutable std::ostringstream info_str_;

  void make_info_str(const input::Parameters& inputs);
  void print_progress(const int& num_measurement, const int& num_measure_steps) const;
};

} // end namespace vmc

#endif

/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:07:37
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-18 12:33:52
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <deque>
#include <tuple>
#include <stdexcept>
#include <Eigen/LU>
#include "../utils/utils.h"
#include "./vmcrun.h"
#include "./opt/prob_linesearch.h"

namespace vmc {

class Optimizer 
{
public:
  Optimizer() {} 
  Optimizer(const input::Parameters& parms); 
  ~Optimizer() {}
  int init(const input::Parameters& parms, const VMCRun& vmc);
  int print_info(const VMCRun& vmc);
  int optimize(VMCRun& vmc);
private:
  enum class exit_status {converged, notconvgd, maxiter, terminated};
  enum class cg_type {FR, PR, DY, HS};
  int num_parms_;
  int num_parms_print_;
  int num_vmc_samples_;
  mcdata::MC_Observable optimal_parms_;
  mcdata::MC_Observable optimal_energy_;
  mcdata::MC_Observable energy_error_bar_;
  var::parm_vector lbound_;
  var::parm_vector ubound_;
  std::vector<double> xvar_values_;
  bool all_converged_{true};

  // MK statistics parameters
  util::MK_Statistic mk_statistic_;
  util::MK_Statistic mk_statistic_en_;
  int mk_series_len_{30};
  double mk_thresold_{0.30};

  // series energy values
  std::deque<double> iter_energy_;
  std::deque<double> iter_energy_err_;
  std::deque<double> iter_gnorm_;

  // Probabilistic line search
  opt::ProbLineSearch probLS_;
  int num_opt_samples_{20};
  int maxiter_{200};
  double start_tstep_{0.1};
  double search_tstep_{0.1};
  double grad_tol_{5.0E-3};
  //double ftol_{5.0E-5};
  bool print_progress_{false};
  bool print_log_{true};

  // Stochastic CG parameters
  cg_type cg_algorithm_{cg_type::PR};
  int cg_maxiter_{0};

  // SR parameters
  bool dir_cutoff_{false};
  double stabilizer_{0.1};
  double w_svd_cut_{0.001};

  // probLS parameters
  int num_probls_steps_{50}; 
  double pls_c1_{0.05};
  double pls_cW_{0.3};
  double pls_alpha0_{0.02}; 
  //double pls_target_df_{0.5};

  // fixed step size 
  int refinement_cycle_{50};
  int refinement_level_{0};
  int fixedstep_iter_{0};


  // progress file
  std::ofstream logfile_;
  std::ofstream file_energy_;
  std::ofstream file_vparms_;
  std::ofstream file_life_;
  std::string life_fname_;
  std::string iter_efile_;
  std::string iter_vfile_;

  int init_sample(const VMCRun& vmc, const int& sample);
  int finalize_sample(const exit_status& status);
  double line_search(VMCRun& vmc, RealVector& vparms, double& en, 
    double& en_err, RealVector& grad, RealVector& grad_err, 
    const RealVector& search_dir);
  int stochastic_reconf(const RealVector& grad, RealMatrix& srmat, RealMatrix& work_mat,
    RealVector& search_dir);
  int stochastic_CG(const RealVector& grad, const RealVector& grad_xprev, 
  RealVector& stochastic_grad, RealVector& search_dir);
  double projected_gnorm(const RealVector& x, const RealVector& grad, 
    const RealVector& lb, const RealVector& ub) const;
  double series_avg(const std::deque<double>& data) const;
  std::ostream& print_progress(std::ostream& os, const int& iter, const std::string& algorithm);
  std::ostream& print_progress(std::ostream& os, 
    const var::parm_vector& vparms, const double& energy,
    const double& error_bar, const RealVector& grad, const RealVector& search_dir,
    const double& gnorm, const double& avg_gnorm, const double& en_trend,
    const double& mk_trend, 
    const int& trend_elem);
};


} // end namespace vmc

#endif

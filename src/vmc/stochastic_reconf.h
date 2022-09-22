/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:07:37
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-18 12:33:52
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef STOCHASTIC_RECONF_H
#define STOCHASTIC_RECONF_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <deque>
#include "../utils/utils.h"
#include "./vmc.h"

namespace vmc {

enum class exit_stat {converged, maxiter, terminated};

class StochasticReconf 
{
public:
  StochasticReconf() {} 
  StochasticReconf(const input::Parameters& parms); 
  ~StochasticReconf() {}
  int init(const input::Parameters& parms, const VMC& vmc);
  int print_info(const VMC& vmc);
  int optimize(VMC& vmc);
  int optimize_old(VMC& vmc);
  const var::parm_vector& optimal_parms(void) const { return vparms_; }
  const int& num_opt_samples(void) const { return num_opt_samples_; }
  int start(const VMC& vmc, const int& sample);
  int iterate(var::parm_vector& vparms, const double& energy,
    const double& energy_err, const Eigen::VectorXd& grad, 
    Eigen::MatrixXd& sr_matrix, exit_stat& status);
  int finalize(void);
private:
  int num_parms_;
  int num_parms_print_;
  mcdata::MC_Observable optimal_parms_;
  mcdata::MC_Observable optimal_energy_;
  mcdata::MC_Observable energy_error_bar_;
  var::parm_vector vparms_;
  var::parm_vector vparms_start_;
  var::parm_vector search_dir_;
  var::parm_vector lbound_;
  var::parm_vector ubound_;
  var::parm_vector range_;
  Eigen::MatrixXd sr_matrix_;
  Eigen::VectorXd grad_;
  std::vector<double> xvar_values_;
  // Mann-Kendall trend test for converegence
  util::MK_Statistic mk_statistic_;

  // series energy values
  std::deque<double> iter_energy_;
  std::deque<double> iter_energy_err_;

  int refinement_cycle_{50};
  int refinement_level_{0};
  int mk_series_len_{30};
  double mk_thresold_{0.30};
  // optimization parameters
  int num_sim_samples_{1000};
  int num_sim_samples2_{5000};
  int num_opt_samples_{20};
  int sample_{0};
  int iter_{0};
  int max_iter_{200};
  int flat_tail_len_{20};
  int random_start_{true};
  bool ln_search_{false};
  bool SA_steps_{true};
  bool dir_cutoff_{false};
  double lambda_cut_{1.0E-3};
  double stabilizer_{0.1};
  double w_svd_cut_{0.001};
  double grad_tol_{5.0E-3};
  double ftol_{5.0E-5};
  double lnsearch_mu_{1.0E-2};
  double lnsearch_beta_{0.5};
  double lnsearch_step_{0.5};
  double start_tstep_{0.1};
  double search_tstep_{0.1};
  //double rslope_tol_{1.0E-2};
  //double aslope_tol_{1.0E-5};
  bool print_progress_{false};
  bool print_log_{true};
  Eigen::VectorXi converged_;


  // progress file
  std::ofstream logfile_;
  std::ofstream file_energy_;
  std::ofstream file_vparms_;
  std::ofstream file_life_;
  std::string life_fname_;
  std::string iter_efile_;
  std::string iter_vfile_;

  std::ostream& print_progress(std::ostream& os, const var::parm_vector& vparms, 
    const double& energy, const double& error_bar, const Eigen::VectorXd& grad, 
    const double& gnorm);
  RealVector lsqfit(const std::vector<double>& iter_energy, const int& max_fitpoints) const;
  void fit_vparms(const std::deque<mcdata::data_t>& iter_vparms, RealVector& fit_slope) const;
  void iter_mean(const std::deque<mcdata::data_t>& iter_vparms, RealVector& fit_slope) const;
  RealVector lsqfit_poly(const std::vector<double>& iter_energy, const int& poly_deg) const;
  bool line_search(VMC& vmc, Eigen::VectorXd& xp, double& fx, double& fx_err,
    Eigen::VectorXd& grad, const Eigen::VectorXd& dir, const double& max_step);
  double max_step_size(const Eigen::VectorXd& x0, const Eigen::VectorXd& dir, 
    const Eigen::VectorXd& lb, const Eigen::VectorXd& ub);
  double proj_grad_norm(const Eigen::VectorXd& x, const Eigen::VectorXd& grad, 
    const Eigen::VectorXd& lb, const Eigen::VectorXd& ub);
  void get_stabilized_dir(const Eigen::MatrixXd& sr_matrix, const Eigen::VectorXd& grad,
    Eigen::VectorXd& dir);
};


} // end namespace vmc

#endif

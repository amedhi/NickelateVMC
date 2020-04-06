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
//#include "../utils/utils.h"
#include "./vmc.h"

namespace vmc {

class StochasticReconf 
{
public:
  StochasticReconf() {} 
  StochasticReconf(const input::Parameters& parms); 
  ~StochasticReconf() {}
  int init(const input::Parameters& parms, const VMC& vmc);
  int optimize(VMC& vmc);
  const var::parm_vector& optimal_parms(void) const { return vparms_; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  int num_parms_;
  mcdata::MC_Observable optimal_parms_;
  mcdata::MC_Observable optimal_energy_;
  mcdata::MC_Observable energy_error_bar_;
  var::parm_vector vparms_;
  var::parm_vector vparms_start_;
  var::parm_vector lbound_;
  var::parm_vector ubound_;
  var::parm_vector range_;
  Eigen::MatrixXd sr_matrix_;
  Eigen::VectorXd grad_;
  std::vector<double> xvar_values_;
  // Mann-Kendall trend test for converegence
  // util::MK_Statistic mk_statistic_;
  // optimization parameters
  int num_sim_samples_{1000};
  int num_sim_samples2_{5000};
  int num_opt_samples_{30};
  int max_iter_{200};
  int flat_tail_len_{20};
  int random_start_{true};
  double stabilizer_{1.0E-4};
  double grad_tol_{0.01};
  double ftol_{0.01};
  //int refinement_cycle_{100};
  //int mk_series_len_{40};
  double lnsearch_mu_{1.0E-4};
  double lnsearch_beta_{0.5};
  double lnsearch_step_{0.5};
  //double start_tstep_{0.05};
  //double rslope_tol_{1.0E-2};
  //double aslope_tol_{1.0E-5};
  //double mk_thresold_{0.30};
  bool print_progress_{false};
  bool print_log_{true};

  // progress file
  std::ofstream logfile_;
  std::ofstream file_energy_;
  std::ofstream file_vparms_;
  std::ofstream file_life_;
  std::string life_fname_;
  std::string iter_efile_;
  std::string iter_vfile_;

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
};


} // end namespace vmc

#endif

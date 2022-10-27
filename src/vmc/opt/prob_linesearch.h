/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:07:37
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-18 12:33:52
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef PROB_LINESEARCH_H
#define PROB_LINESEARCH_H

#include <iostream>
#include <fstream>
#include <stdexcept>
#include "../../wavefunction/matrix.h"
#include <Eigen/LU>

namespace opt {

class ProbLineSearchGP
{
public:
  ProbLineSearchGP(); 
  ~ProbLineSearchGP() {}
  int set_parameters(const double& theta, const double& offset);
  int add_obs(const double& t, const double& f, const double& df,
    const double& fvar=0.0, const double& dfvar=0.0);
  int reset(void);
  int update(void);
  double expected_improvement(const double& t);
  std::vector<double> quad_poly_solve(const double& a, const double& b, 
    const double& c, const double& val) const;
  int test(void);
private:
  bool ready_{false};
  // Hyperparamters of the GP
  double theta_{1.0};
  double theta_sq_{1.0};
  double offset_{10.0};
  double PI_{3.141592653589793};

  // Observation counter and arrays to store observations 
  int N_{0};
  bool have_min_obs_{false};
  double min_posterior_mean_;
  std::vector<double> ts_;
  std::vector<double> fs_;
  std::vector<double> dfs_;
  std::vector<double> fvars_;
  std::vector<double> dfvars_;

  // Kernel matrices
  RealMatrix kernel_K_;
  RealMatrix kernel_Kd_;
  RealMatrix kernel_dKd_;

  // Gram matrix and pre-computed "weights" of the GP
  RealMatrix G_;
  RealVector w_;

  //mutable Eigen::FullPivLU<RealMatrix> FullPivLU_; 

public:
  RealVector solve_G(const RealVector& b) const;
  double mu(const double& t) const;
  double dmu(const double& t) const;
  double d2mu(const double& t) const;
  double d3mu(const double& t) const;
  double V(const double& t) const;
  double Vd(const double& t) const;
  double dVd(const double& t) const;
  double Cov_0(const double& t) const;
  double Covd_0(const double& t) const;
  double dCov_0(const double& t) const;
  double dCovd_0(const double& t) const;
  int cubic_poly_coeff(const double& t, double& a, double& b, double& c, double& d) const;
  int quad_poly_coeff(const double& t, double& a, double& b, double& c) const;
  std::vector<double> find_dmu_equal(const double& val) const;
  std::vector<double> find_cubic_minima(void) const;
  double get_k(const double& x, const double& y) const;
  double get_kd(const double& x, const double& y) const;
  double get_dkd(const double& x, const double& y) const;
  double get_d2k(const double& x, const double& y) const;
  double get_d3k(const double& x, const double& y) const;
  double get_d2kd(const double& x, const double& y) const;
};

class ProbLineSearch
{
public:
  using fquad = std::tuple<double,double,double,double>;
  ProbLineSearch();
  ~ProbLineSearch() {}
  void set_parameters(const double& c1=0.05, const double& cW=0.3, 
    const double& fpush=1.0, const double& alpha0=0.01, 
    const double& target_df=0.5, const double& df_lo=-0.1, 
    const double& df_hi=1.1, const int& max_steps=10, 
    const int& max_expl=6, const double& max_dmu0=0.0, 
    const double& max_change_factor=10.0, 
    const std::string& expl_policy="linear");
  void set_c1(const double& c1) { c1_=c1; }
  void set_cW(const double& cW) { cW_=cW; }
  void set_alpha0(const double& alpha0) { alpha0_=alpha0; alpha_stats_=alpha0; }
  void set_target_df(const double& df) { target_df_=df; }
  void set_fpush(const double& fpush) { fpush_=fpush; }
  double start(const double& en, const double& en_err, const RealVector& grad, 
    const RealVector& grad_err, const RealVector& search_dir);
  double do_step(const double& en, const double& en_err, const RealVector& grad, const RealVector& grad_err, 
    const RealVector& search_dir, bool& accept_status, std::string& message);
  int test(void);
private:
  ProbLineSearchGP GP_;
  double c1_;
  double cW_; 
  double fpush_; 
  double alpha0_; 
  double alpha_stats_; 
  double target_df_;
  double df_lo_;
  double df_hi_; 
  int max_steps_; 
  int max_expl_;
  double max_dmu0_; 
  double max_change_factor_;
  std::string expl_policy_;

  // projected quantities
  double f0_;
  double df0_;

  // running parameters
  bool new_search_{false};
  int num_steps_{0};
  int num_expl_{0};
  double current_t_{1.0};
  bool abort_status1_{false};
  bool abort_status2_{false};


double variance(const double& ferr);
double projected_grad(const RealVector& grad, const RealVector& search_dir);
double projected_gradvar(const RealVector& grad_err, const RealVector& search_dir);
fquad scale_obs(const double& raw_f, const double& raw_df, 
  const double& raw_fvar, const double& raw_dfvar);
double rescale_t(const double& raw_t) const;
double find_next_t(void);
double find_abort_t(void);
double compute_p_wolfe(const double& t);
bool check_acceptibility(void);
double bounded_bivariate_normal_integral(const double& rho, const double& xl, 
  const double& xu, const double& yl, const double& yu) const;
double unbounded_bivariate_normal_integral(const double& rho, const double& xl, 
  const double& xu) const;
double CDF(const double& x) const;
};


} // end namespace vmc

#endif

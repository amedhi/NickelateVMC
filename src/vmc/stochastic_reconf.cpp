/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:19:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 17:40:34
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "./stochastic_reconf.h"

namespace vmc {

int StochasticReconf::init(const input::Parameters& inputs, const VMC& vmc) 
{
  // problem size
  num_parms_ = vmc.num_varp();
  vparms_.resize(num_parms_);
  vparms_start_.resize(num_parms_);
  lbound_ = vmc.varp_lbound();
  ubound_ = vmc.varp_ubound();
  range_ = ubound_-lbound_;
  grad_.resize(num_parms_);
  sr_matrix_.resize(num_parms_,num_parms_);

  // optimization parameters
  int nowarn;
  num_sim_samples_ = inputs.set_value("measure_steps", -1, nowarn);
  num_sim_samples2_ = 5*num_sim_samples_;
  num_opt_samples_ = inputs.set_value("sr_opt_samples", 10, nowarn);
  max_iter_ = inputs.set_value("sr_max_iter", 100, nowarn);
  //start_tstep_ = inputs.set_value("sr_start_tstep", 0.05, nowarn);
  random_start_ = inputs.set_value("sr_random_start", true, nowarn);
  flat_tail_len_ = inputs.set_value("sr_flat_tail_len", 10, nowarn);
  ftol_ = inputs.set_value("sr_ftol", 5.0E-5, nowarn);
  grad_tol_ = inputs.set_value("sr_gradtol", 5.0E-3, nowarn);
  print_progress_ = inputs.set_value("sr_progress_stdout", false, nowarn);
  print_log_ = inputs.set_value("sr_progress_log", true, nowarn);

  //rslope_tol_ = inputs.set_value("sr_rslope_tol", 1.0E-2, nowarn);
  //mk_series_len_ = inputs.set_value("sr_series_len", 40, nowarn);
  //stabilizer_ = inputs.set_value("sr_stabilizer", 1.0E-4, nowarn);
  //mk_thresold_ = inputs.set_value("sr_fluctuation_tol", 0.30, nowarn);

  // optimal parameter values
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  bool replace_mode = true;
  if (mode=="APPEND") replace_mode = false;
  optimal_parms_.init("opt_params");
  optimal_parms_.set_ofstream(vmc.prefix_dir());
  optimal_parms_.resize(vmc.num_varp(), vmc.varp_names());
  optimal_parms_.set_replace_mode(replace_mode);
  optimal_parms_.switch_on();
  optimal_parms_.switch_on();
  std::vector<std::string> elem_names{"energy"};
  optimal_energy_.init("opt_energy");
  optimal_energy_.set_ofstream(vmc.prefix_dir());
  optimal_energy_.resize(1, elem_names);
  optimal_energy_.set_replace_mode(replace_mode);
  optimal_energy_.switch_on();
  energy_error_bar_.init("error_bar");
  energy_error_bar_.switch_on();

  // observable file header
  std::stringstream heading;
  vmc.print_info(heading);
  optimal_parms_.print_heading(heading.rdbuf()->str(), vmc.xvar_names());
  optimal_energy_.print_heading(heading.rdbuf()->str(), vmc.xvar_names());
  xvar_values_ = vmc.xvar_values();

  // iteration related file names
  life_fname_ = vmc.prefix_dir()+"ALIVE.d";
  iter_efile_ = vmc.prefix_dir()+"iter_energy.txt";
  iter_vfile_ = vmc.prefix_dir()+"iter_params.txt";

  // progress file
  if (print_log_) {
    logfile_.open(vmc.prefix_dir()+"log_optimization.txt");
    if (!logfile_.is_open())
      throw std::runtime_error("StochasticReconf::init: file open failed");
    vmc.print_info(logfile_);
    logfile_ << "#" << std::string(72, '-') << std::endl;
    logfile_ << "Stochastic Reconfiguration" << std::endl;
    logfile_ << "num_opt_samples = " << num_opt_samples_ << std::endl;
    logfile_ << "max_iter = " << max_iter_ << std::endl;
    logfile_ << "random_start = " << random_start_ << std::endl;
    logfile_ << "stabilizer = " << stabilizer_ << std::endl;
    logfile_ << "flat_tail_len = " << flat_tail_len_ << std::endl;
    logfile_ << "grad_tol = " << grad_tol_ << std::endl;
    logfile_ << "ftol = " << ftol_ << std::endl;
    logfile_ << "#" << std::string(72, '-') << std::endl;
  }
  if (print_progress_) {
    std::cout << "#" << std::string(72, '-') << std::endl;
    std::cout << "Stochastic Reconfiguration" << std::endl;
    std::cout << "num_opt_samples = " << num_opt_samples_ << std::endl;
    std::cout << "max_iter = " << max_iter_ << std::endl;
    std::cout << "random_start = " << random_start_ << std::endl;
    std::cout << "stabilizer = " << stabilizer_ << std::endl;
    std::cout << "flat_tail_len = " << flat_tail_len_<< std::endl;
    std::cout << "grad_tol = " << grad_tol_ << std::endl;
    std::cout << "ftol = " << ftol_ << std::endl;
    std::cout << "#" << std::string(72, '-') << std::endl;
  }
  return 0;
}

int StochasticReconf::optimize(VMC& vmc)
{
  /*
  //std::vector<double> en = {0.3736,2.4375,3.9384,3.3123,5.4947,5.4332,6.3932,9.0605,9.3642,9.5207};
  std::vector<double> en = {-0.0546347,-0.0555617,-0.0559485,-0.0555145,-0.0560293,-0.0566849,
    -0.0571339,-0.0576311,-0.0584373,-0.0582745,-0.0586636,-0.0585779,
    -0.0587886,-0.0589912,-0.0589437,-0.0596863,-0.0593831,-0.0597698,
    -0.0593647,-0.0593304,-0.0603898,-0.0599514,-0.0592811,-0.0596215};
  if_converged(en);
  */
  xvar_values_ = vmc.xvar_values();
  optimal_parms_.reset();
  optimal_energy_.reset();
  energy_error_bar_.reset();
  // start optimization
  // starting value of variational parameters
  vparms_start_ = vmc.varp_values();
  bool all_samples_converged = true;
  for (int sample=1; sample<=num_opt_samples_; ++sample) {
    // iteration files
    std::ofstream life_fs(life_fname_);
    life_fs.close();
    file_energy_.open(iter_efile_);
    file_vparms_.open(iter_vfile_);
    vmc.print_info(file_energy_);
    file_energy_<<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
    file_vparms_<<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
    file_energy_ << "# Results: " << "Iteration Energy";
    file_energy_ << " (sample "<<sample<<" of "<<num_opt_samples_<<")\n";
    file_energy_ << "#" << std::string(72, '-') << "\n";
    file_energy_ << "# ";
    file_energy_ << std::left;
    file_energy_ <<std::setw(7)<<"iter"<<std::setw(15)<<"energy"<<std::setw(15)<<"err";
    file_energy_ << std::endl;
    file_energy_ << "#" << std::string(72, '-') << "\n";
    file_energy_ << std::right;

    vmc.print_info(file_vparms_);
    file_vparms_ << "# Results: " << "Iteration Variational Parameters";
    file_vparms_ << " (sample "<<sample<<" of "<<num_opt_samples_<<")\n";
    file_vparms_ << "#" << std::string(72, '-') << "\n";
    file_vparms_ << "# ";
    file_vparms_ << std::left;
    file_vparms_ << std::setw(7)<<"iter";
    for (const auto& p : vmc.varp_names()) file_vparms_<<std::setw(15)<<p.substr(0,14);
    file_vparms_ << std::endl;
    file_vparms_ << "#" << std::string(72, '-') << "\n";
    file_vparms_ << std::right;

    // message
    if (print_progress_) {
      std::cout << " Starting sample "<<sample<<" of "<<num_opt_samples_<<"\n"; 
      std::cout << std::flush;
    }
    if (print_log_) {
      logfile_ << " Starting sample "<<sample<<" of "<<num_opt_samples_<<"\n"; 
      logfile_ << std::flush;
    }

    // give gaussian randomness to the parameters
    if (random_start_) {
      for (int i=0; i<num_parms_; ++i) {
        double sigma = 0.10*(ubound_[i]-lbound_[i]);
        std::normal_distribution<double> random_normal(vparms_start_[i],sigma);
        double x = random_normal(vmc.rng());
        if (x>=lbound_[i] && x<=ubound_[i]) {
          vparms_[i] = x;
        }
        else {
          vparms_[i] = vparms_start_[i];
        }
      }
    }
    else {
      vparms_ = vparms_start_;
    }
    // seed vmc rng
    int rng_seed = 1;
    //vmc.seed_rng(1);
    // initial function values
    double start_en, error_bar; 
    vmc.sr_function(vparms_,start_en,error_bar,grad_,sr_matrix_,num_sim_samples_,rng_seed);
    // apply to stabilizer to sr matrix 
    for (int i=0; i<num_parms_; ++i) sr_matrix_(i,i) += stabilizer_;

    // Stochastic reconfiguration iterations
    bool converged = false;
    double energy = start_en;
    std::vector<double> iter_energy;
    std::vector<double> iter_energy_err;
    std::vector<double> iter_gnorm;
    std::deque<mcdata::data_t> iter_vparms;
    std::deque<mcdata::data_t> iter_grads;
    RealVector vparm_slope(num_parms_);
    RealVector grads_slope(num_parms_);
    bool ln_search = true;
    for (int iter=1; iter<=max_iter_; ++iter) {
      // file outs
      file_energy_<<std::setw(6)<<iter<<std::scientific<<std::setw(16)<<energy; 
      file_vparms_<<std::setw(6)<<iter; 
      // progress log
      if (print_progress_) {
        std::ios state(NULL);
        state.copyfmt(std::cout);
        std::cout << "#" << std::string(60, '-') << std::endl;
        std::cout << std::left;
        if (ln_search) 
          std::cout << " iteration  =  "<<std::setw(6)<<iter<<"(line search)"<<"\n";
        else
          std::cout << " iteration  =  "<<std::setw(6)<<iter<<"(SA step)"<<"\n";
      }
      if (print_log_) {
        logfile_ << "#" << std::string(60, '-') << std::endl;
        logfile_ << std::left;
        if (ln_search) 
          logfile_ << " iteration  =  "<<std::setw(6)<<iter<<"(line search)"<<"\n";
        else 
          logfile_ << " iteration  =  "<<std::setw(6)<<iter<<"(SA step)"<<"\n";
      }

      // search direction
      Eigen::VectorXd search_dir = sr_matrix_.fullPivLu().solve(-grad_);

      // max step size considering bounds
      double max_step = max_step_size(vparms_, search_dir, lbound_, ubound_);
      if (ln_search) {
        // line search with Armijio condition
        ln_search = line_search(vmc,vparms_,energy,error_bar,grad_,search_dir,max_step);
        // simulation at new parameters
        vmc.sr_function(vparms_,energy,error_bar,grad_,sr_matrix_,num_sim_samples_,rng_seed);
      }
      else {
        // Stochastic Approximation steps
        double tstep = 1.0/(iter+1);
        vparms_ += tstep * search_dir;
        // box constraint and max_norm (of components not hitting boundary) 
        vparms_ = lbound_.cwiseMax(vparms_.cwiseMin(ubound_));
        // simulation at new parameters with higher precision
        vmc.sr_function(vparms_,energy,error_bar,grad_,sr_matrix_,num_sim_samples2_,rng_seed);
        // apply to stabilizer to sr matrix 
        for (int i=0; i<num_parms_; ++i) sr_matrix_(i,i) += stabilizer_;
      }
      // simulation at new parameters
      //vmc.sr_function(vparms_,en,error_bar,grad_,sr_matrix_,num_sim_samples_,rng_seed);
      double gnorm = grad_.squaredNorm();
      double proj_norm = proj_grad_norm(vparms_, grad_, lbound_, ubound_);
      gnorm = std::min(gnorm, proj_norm);

      // file outs
      file_energy_<<std::fixed<<std::setw(10)<<error_bar<<std::endl<<std::flush;
      for (int i=0; i<vparms_.size(); ++i) file_vparms_<<std::setw(15)<<vparms_[i];
      file_vparms_<<std::endl<<std::flush;
      // progress log
      if (print_progress_) {
        std::ios state(NULL);
        state.copyfmt(std::cout);
        std::cout <<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
        std::cout << " grad       =";
        for (int i=0; i<grad_.size(); ++i)
          std::cout << std::setw(15)<<grad_[i];
        std::cout << "\n";
        std::cout << " search_dir =";
        for (int i=0; i<search_dir.size(); ++i)
          std::cout << std::setw(15)<<search_dir[i];
        std::cout << "\n";
        std::cout << " varp       =";
        for (int i=0; i<vparms_.size(); ++i)
          std::cout << std::setw(15)<<vparms_[i];
        std::cout << "\n";
        std::cout << " gnorm      ="<<std::setw(15)<< gnorm<<"\n";
        std::cout << " energy     ="<<std::setw(15)<<energy << "   (+/-) ";
        std::cout <<std::fixed<<std::setw(10)<<error_bar<<"\n";
        std::cout.copyfmt(state);
      }
      if (print_log_) {
        logfile_ <<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
        logfile_ << " grad       =";
        for (int i=0; i<grad_.size(); ++i)
          logfile_ << std::setw(15)<<grad_[i];
        logfile_ << "\n";
        logfile_ << " search_dir =";
        for (int i=0; i<search_dir.size(); ++i)
          logfile_ << std::setw(15)<<search_dir[i];
        logfile_ << "\n";
        logfile_ << " varp       =";
        for (int i=0; i<vparms_.size(); ++i)
          logfile_ << std::setw(15)<<vparms_[i];
        logfile_ << "\n";
        logfile_ << " gnorm      ="<<std::setw(15)<< gnorm<<"\n";
        logfile_ << " energy     ="<<std::setw(15)<<energy << "   (+/-) "; 
        logfile_ << std::fixed<<std::setw(10)<<error_bar<<"\n";
      }

      // check convergence
      if (!ln_search) {
        iter_energy.push_back(energy);
        iter_energy_err.push_back(error_bar);
        iter_gnorm.push_back(gnorm);
        iter_vparms.push_back(vparms_);
        int len = iter_energy.size();
        if (len > flat_tail_len_) {
          iter_vparms.pop_front();
          RealVector efit = lsqfit(iter_energy,flat_tail_len_);
          double slope = efit[1];
          double gnorm_avg = 0.0;
          for (int i=1; i<=flat_tail_len_; ++i) gnorm_avg += iter_gnorm[len-i];
          gnorm_avg /= flat_tail_len_;
          //std::cout << "gnorm_avg = " << gnorm_avg << "\n";
          //std::cout << "fslope    = " << slope << "\n";
          if (gnorm_avg<grad_tol_ && slope<ftol_) {
            converged = true; 
            break;
          }
          else {
            converged = false;
          }
        }
      }
      if (!boost::filesystem::exists(life_fname_)) break;
    } // iteration over
    if (!converged) all_samples_converged = false;
    if (print_progress_) {
      std::cout << "#" << std::string(60, '-') << std::endl;
      if (converged) std::cout<<" Converged!"<<std::endl;
      else std::cout <<" NOT converged"<<std::endl;
      std::cout << "#" << std::string(60, '-') << std::endl << std::flush;
    }
    if (print_log_) {
      logfile_ << "#" << std::string(60, '-') << std::endl;
      if (converged) logfile_<<" Converged!"<<std::endl;
      else logfile_ <<" NOT converged"<<std::endl;
      logfile_ << "#" << std::string(60, '-')<<std::endl<<std::endl<< std::flush;
    }
    file_energy_.close();
    file_vparms_.close();

    // optimized quantities for this sample
    if (converged) {
      int n = iter_energy.size()-flat_tail_len_;
      for (int i=0; i<flat_tail_len_; ++i) {
        optimal_energy_ << iter_energy[n+i];
        energy_error_bar_ << iter_energy_err[n+i];
        optimal_parms_ << iter_vparms[i];
      }
      // print sample values
      // optimal energy
      optimal_energy_.open_file();
      optimal_energy_.fs()<<std::right;
      optimal_energy_.fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
      optimal_energy_.fs()<<"#";
      for (const auto& p : xvar_values_) 
        optimal_energy_.fs()<<std::setw(13)<<p;
      optimal_energy_.fs()<<std::setw(15)<<optimal_energy_.mean();
      optimal_energy_.fs()<<std::fixed<<std::setw(10)<<energy_error_bar_.mean();
      optimal_energy_.fs()<<std::setw(10)<<num_sim_samples2_; 
      if (converged) 
        optimal_energy_.fs()<<std::setw(11)<<"CONVERGED"<<std::setw(7)<<0<<std::endl;
      else 
        optimal_energy_.fs()<<std::setw(11)<<"NOT_CONVD"<<std::setw(7)<<0<<std::endl;
      // save the values
      optimal_energy_.save_result();
      energy_error_bar_.save_result();
      // grand average for samples done so far
      optimal_energy_.fs()<<std::right;
      optimal_energy_.fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
      optimal_energy_.fs() << "#" << std::string(72, '-') << "\n";
      optimal_energy_.fs() << "# grand average:\n"; // 
      optimal_energy_.fs() << "#" << std::string(72, '-') << "\n";
      for (const auto& p : xvar_values_) 
        optimal_energy_.fs()<<std::setw(14)<<p;
      optimal_energy_.fs()<<std::setw(15)<<optimal_energy_.grand_data().mean();
      optimal_energy_.fs()<<std::fixed<<std::setw(10)<<energy_error_bar_.grand_data().mean();
      optimal_energy_.fs()<<std::setw(10)<<num_sim_samples2_; 
      if (all_samples_converged) 
        optimal_energy_.fs()<<std::setw(11)<<"CONVERGED"<<std::setw(7)<<0<<std::endl;
      else 
        optimal_energy_.fs()<<std::setw(11)<<"NOT_CONVD"<<std::setw(7)<<0<<std::endl;
      optimal_energy_.close_file();
      //optimal_energy_.print_result(xvar_values_);
      // optimal variational parameters
      optimal_parms_.open_file();
      optimal_parms_.fs()<<"#";
      optimal_parms_.print_result(xvar_values_);
      // save the values
      optimal_parms_.save_result();
      // grand average for samples done so far
      optimal_parms_.print_grand_result(xvar_values_);

      // reset for next sample
      optimal_energy_.reset();
      energy_error_bar_.reset();
      optimal_parms_.reset();
    } // one sample converged
  } // samples
  logfile_.close();
  // print results

  // grand averages
  // optimal energy
  optimal_energy_.open_file();
  optimal_energy_.fs()<<std::right;
  optimal_energy_.fs()<<std::scientific<<std::uppercase<<std::setprecision(6);
  optimal_energy_.fs() << "#" << std::string(72, '-') << "\n";
  optimal_energy_.fs() << "# grand average:\n"; // 
  optimal_energy_.fs() << "#" << std::string(72, '-') << "\n";
  for (const auto& p : xvar_values_) 
    optimal_energy_.fs()<<std::setw(14)<<p;
  optimal_energy_.fs()<<std::setw(15)<<optimal_energy_.grand_data().mean();
  optimal_energy_.fs()<<std::fixed<<std::setw(10)<<energy_error_bar_.grand_data().mean();
  optimal_energy_.fs()<<std::setw(10)<<num_sim_samples_; 
  if (all_samples_converged) 
    optimal_energy_.fs()<<std::setw(11)<<"CONVERGED"<<std::setw(7)<<0<<std::endl;
  else 
    optimal_energy_.fs()<<std::setw(11)<<"NOT_CONVD"<<std::setw(7)<<0<<std::endl;
  optimal_energy_.close_file();
  optimal_parms_.print_grand_result(xvar_values_);
  // reset
  optimal_energy_.reset_grand_data();
  energy_error_bar_.reset_grand_data();
  optimal_parms_.reset_grand_data();

  return 0;
}

RealVector StochasticReconf::lsqfit(const std::vector<double>& iter_energy, 
  const int& max_fitpoints) const
{
  /* 
     Least Square Fit with 'f(x) = p0 + p1*x' (by QR method) 
  */
  int num_points = max_fitpoints;
  if (num_points > iter_energy.size()) num_points = iter_energy.size();

  //  The 'xrange' is set at [0:1.0]
  double dx = 1.0/num_points;
  RealMatrix A(num_points,2);
  for (int i=0; i<num_points; ++i) {
    A(i,0) = 1.0;
    A(i,1) = (i+1) * dx;
  }
  // take the 'last num_points' values
  int n = iter_energy.size()-num_points;
  RealVector b(num_points);
  for (int i=0; i<num_points; ++i) b[i] = iter_energy[n+i];
  RealVector p = A.colPivHouseholderQr().solve(b);
  //std::cout << "poly = " << p.transpose() << "\n"; 
  //std::cout << "slope = " << p(1) << "\n"; 
  //getchar();
  return p;
}

void StochasticReconf::fit_vparms(const std::deque<mcdata::data_t>& iter_vparms, RealVector& fit_slope) const
{
  int num_points = iter_vparms.size();
  //  The 'xrange' is set at [0:1.0]
  double dx = 1.0/num_points;
  RealMatrix A(num_points,2);
  for (int i=0; i<num_points; ++i) {
    A(i,0) = 1.0;
    A(i,1) = (i+1) * dx;
  }
  for (int n=0; n<num_parms_; ++n) {
    RealVector b(num_points);
    for (int i=0; i<num_points; ++i) b[i] = iter_vparms[i][n];
    RealVector p = A.colPivHouseholderQr().solve(b);
    fit_slope(n) = p(1);
  }
} 

void StochasticReconf::iter_mean(const std::deque<mcdata::data_t>& iter_vparms, RealVector& mean) const
{
  int num_points = iter_vparms.size();
  mean.setZero();
  for (auto& elem : iter_vparms) {
    for (int i=0; i<num_parms_; ++i) {
      mean(i) += elem(i);
    }
  }
  for (int i=0; i<num_parms_; ++i) {
    mean(i) /= num_points;
  }
} 

RealVector StochasticReconf::lsqfit_poly(const std::vector<double>& iter_energy, 
  const int& poly_deg) const
{
  /* 
     Least Square Fit (by QR method) of the data with 
     polynomial of degree 'fitpoly_degree_'
  */
  int num_points = iter_energy.size();
  int poly_degree = 3;
  //  The 'xrange' is set at [0:1.0]
  double dx = 1.0/num_points;
  RealMatrix A(num_points,2);
  for (int i=0; i<num_points; ++i) {
    A(i,0) = 1.0;
    double x = (i+1) * dx;
    for (int n=1; n<=poly_degree; ++n) A(i,n) = A(i,n-1)*x;
  }
  RealVector b(num_points);
  for (int i=0; i<num_points; ++i) b[i] = iter_energy[i];
  RealVector p = A.colPivHouseholderQr().solve(b);
  return p;
}

bool StochasticReconf::line_search(VMC& vmc, Eigen::VectorXd& x,
  double& fx, double& fx_err, Eigen::VectorXd& grad, 
  const Eigen::VectorXd& dir, const double& max_step) 
{
  // save the function value at the current x
  const Eigen::VectorXd x_init = x;
  const double fx_init = fx;
  // projection of gradient on the search direction
  const double dg_init = grad.dot(dir);
  if (dg_init > 0)
     throw std::logic_error("StochasticReconf::line_search: the moving direction uphill\n");

  double decr_factor = lnsearch_mu_*dg_init; 
  double step = 0.8*std::min(lnsearch_step_, max_step);
  int max_linesearch = 5;
  for(int iter = 0; iter < max_linesearch; ++iter) {
    x.noalias() = x_init + step * dir;
    // evaluate this candidate
    vmc(x,fx,fx_err,grad);
    if (fx < fx_init + step*decr_factor) return true;
    step *= lnsearch_beta_;
  }
  // failed: restore initial values
  x.noalias() = x_init;
  fx = fx_init;
  return false;
}

double StochasticReconf::max_step_size(const Eigen::VectorXd& x0, const Eigen::VectorXd& dir, 
  const Eigen::VectorXd& lb, const Eigen::VectorXd& ub)
{
  const int n = x0.size();
  double step = std::numeric_limits<double>::infinity();
  for (int i=0; i<n; ++i) {
    if(dir[i] > 0.0) {
      step = std::min(step, (ub[i] - x0[i])/dir[i]);
    } 
    else if(dir[i] < 0.0) {
      step = std::min(step, (lb[i] - x0[i])/dir[i]);
    }
  }
  return step;
}

double StochasticReconf::proj_grad_norm(const Eigen::VectorXd& x, const Eigen::VectorXd& grad, 
  const Eigen::VectorXd& lb, const Eigen::VectorXd& ub)
{
  const int n = x.size();
  double res = 0.0;
  for (int i=0; i<n; ++i) {
    double proj = std::max(lb[i], x[i] - grad[i]);
    proj = std::min(ub[i], proj);
    res = std::max(res, std::abs(proj - x[i]));
  }
  return res;
}


} // end namespace vmc

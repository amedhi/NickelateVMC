/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:19:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 17:40:34
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include <limits>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "./optimizer.h"

namespace vmc {

int Optimizer::init(const input::Parameters& inputs, const VMCRun& vmc) 
{
  // Number of variational parameters & bounds
  num_parms_ = vmc.num_varp();
  num_parms_print_ = std::min(num_parms_,10);
  lbound_ = vmc.varp_lbound();
  ubound_ = vmc.varp_ubound();

  //grad_.resize(num_parms_);
  //sr_matrix_.resize(num_parms_,num_parms_);
  //sr_matrix_inv_.resize(num_parms_,num_parms_);

  // optimization parameters
  int nowarn;
  num_opt_samples_ = inputs.set_value("opt_num_samples", 10, nowarn);
  maxiter_ = inputs.set_value("opt_maxiter", 100, nowarn);
  cg_maxiter_ = inputs.set_value("opt_cg_maxiter", 0, nowarn);
  if (cg_maxiter_>0) {
    std::string str = inputs.set_value("opt_cg_algorithm", "PR", nowarn);
    boost::to_upper(str);
    if (str=="FR") {
      cg_algorithm_ = cg_type::FR;
    }
    else if (str=="PR") {
      cg_algorithm_ = cg_type::PR;
    }
    else if (str=="DY") {
      cg_algorithm_ = cg_type::DY;
    }
    else if (str=="HS") {
      cg_algorithm_ = cg_type::HS;
    }
    else {
      throw std::range_error("Optimizer::init: invalid CG method");
    }
  }
  start_tstep_ = inputs.set_value("opt_start_tstep", 0.1, nowarn);
  refinement_cycle_ = inputs.set_value("opt_refinement_cycle", 40, nowarn);
  stabilizer_ = inputs.set_value("opt_stabilizer", 0.2, nowarn);
  w_svd_cut_ = inputs.set_value("opt_svd_cut", 0.001, nowarn);
  //random_start_ = inputs.set_value("sr_random_start", false, nowarn);
  grad_tol_ = inputs.set_value("opt_gradtol", 5.0E-3, nowarn);
  print_progress_ = inputs.set_value("opt_progress_stdout", false, nowarn);
  print_log_ = inputs.set_value("opt_progress_log", true, nowarn);

  // probLS parameters
  num_probls_steps_ = inputs.set_value("opt_pls_steps", 50, nowarn);
  pls_c1_ = inputs.set_value("opt_pls_c1", 0.05, nowarn);
  pls_cW_ = inputs.set_value("opt_pls_cW", 0.3, nowarn);
  pls_alpha0_ = inputs.set_value("opt_pls_alpha0", 0.02, nowarn);
  //pls_target_df_ = inputs.set_value("opt_pls_target_df", 0.5, nowarn);
  probLS_.set_c1(pls_c1_);
  probLS_.set_cW(pls_cW_);
  probLS_.set_alpha0(pls_alpha0_);
  //probLS_.set_target_df(pls_target_df_);

  // Mann-Kendall statistic
  mk_series_len_ = inputs.set_value("opt_mkseries_len", 20, nowarn);
  mk_thresold_ = inputs.set_value("opt_mkfluct_tol", 0.30, nowarn);
  mk_statistic_.resize(num_parms_);
  mk_statistic_en_.resize(1); // for energy
  mk_statistic_.set_maxlen(mk_series_len_);
  mk_statistic_en_.set_maxlen(std::min(6,mk_series_len_));

  // Observables for optimal parameters
  num_vmc_samples_ = vmc.num_measure_steps();
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

  // iteration related files 
  boost::filesystem::path iter_prefix(vmc.prefix_dir()+"/iterations");
  boost::filesystem::create_directory(iter_prefix);
  life_fname_ = vmc.prefix_dir()+"/iterations/ALIVE.d";
  iter_efile_ = vmc.prefix_dir()+"/iterations/iter_energy";
  iter_vfile_ = vmc.prefix_dir()+"/iterations/iter_params";

  return 0;
}

int Optimizer::print_info(const VMCRun& vmc)
{
  if (print_log_) {
    logfile_.open(vmc.prefix_dir()+"log_optimization.txt");
    if (!logfile_.is_open())
      throw std::runtime_error("Optimizer::init: file open failed");
    vmc.print_info(logfile_);
    logfile_ << "#" << std::string(72, '-') << std::endl;
    logfile_ << "Stochastic Reconfiguration" << std::endl;
    logfile_ << "num_opt_samples = " << num_opt_samples_ << std::endl;
    logfile_ << "max_iter = " << maxiter_ << std::endl;
    //logfile_ << "random_start = " << random_start_ << std::endl;
    logfile_ << "stabilizer = " << stabilizer_ << std::endl;
    logfile_ << "start_tstep = " << start_tstep_ << std::endl;
    logfile_ << "w_svd_cut = " << w_svd_cut_ << std::endl;
    logfile_ << "grad_tol = " << grad_tol_ << std::endl;
    //logfile_ << "ftol = " << ftol_ << std::endl;
    logfile_ << "#" << std::string(72, '-') << std::endl;
  }
  if (print_progress_) {
    std::cout << "#" << std::string(72, '-') << std::endl;
    std::cout << "Stochastic Reconfiguration" << std::endl;
    std::cout << "num_opt_samples = " << num_opt_samples_ << std::endl;
    std::cout << "max_iter = " << maxiter_ << std::endl;
    //std::cout << "random_start = " << random_start_ << std::endl;
    std::cout << "stabilizer = " << stabilizer_ << std::endl;
    std::cout << "start_tstep = " << start_tstep_ << std::endl;
    std::cout << "w_svd_cut = " << w_svd_cut_ << std::endl;
    std::cout << "grad_tol = " << grad_tol_ << std::endl;
    //std::cout << "ftol = " << ftol_ << std::endl;
    std::cout << "#" << std::string(72, '-') << std::endl;
  }
  return 0;
}

int Optimizer::init_sample(const VMCRun& vmc, const int& sample)
{
  // reset parameters
  search_tstep_ = start_tstep_;
  fixedstep_iter_ = 1;
  refinement_level_ = 0;

  // Mann-Kendall statistic for converegence test
  mk_statistic_.reset();
  mk_statistic_en_.reset();
  iter_energy_.clear();
  iter_energy_err_.clear();
  iter_gnorm_.clear();

  // iteration files
  std::ostringstream suff;
  suff << "_sample"<<std::setfill('0')<<std::setw(2)<<sample<<".txt";
  file_life_.open(life_fname_);
  file_life_.close();
  file_energy_.open(iter_efile_+suff.str());
  file_vparms_.open(iter_vfile_+suff.str());

  // init Observables for optimal energy & parameters
  xvar_values_ = vmc.xvar_values();
  optimal_parms_.reset();
  optimal_energy_.reset();
  energy_error_bar_.reset();

  vmc.print_info(file_energy_);
  file_energy_<<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
  file_vparms_<<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
  file_energy_ << "# Results: " << "Iteration Energy";
  file_energy_ << " (sample "<<sample<<" of "<<num_opt_samples_<<")\n";
  file_energy_ << "#" << std::string(72, '-') << "\n";
  file_energy_ << "# ";
  file_energy_ << std::left;
  file_energy_ <<std::setw(7)<<"iter"<<std::setw(15)<<"energy"<<std::setw(15)<<"err";
  file_energy_ << "\n";
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
  file_vparms_ << "\n";
  file_vparms_ << "#" << std::string(72, '-') << "\n";
  file_vparms_ << std::right;
  // message
  if (print_progress_) {
    std::cout << " Starting sample "<<sample<<" of "<<num_opt_samples_<<"\n"; 
    //std::cout << "#" << std::string(60, '-') << std::endl;
    std::cout << std::flush;
  }
  if (print_log_) {
    logfile_ << " Starting sample "<<sample<<" of "<<num_opt_samples_<<"\n"; 
    //logfile_ << "#" << std::string(60, '-') << std::endl;
    logfile_ << std::flush;
  }
  return 0;
}

int Optimizer::optimize(VMCRun& vmc)
{
  // Observables
  xvar_values_ = vmc.xvar_values();
  optimal_parms_.reset();
  optimal_energy_.reset();
  energy_error_bar_.reset();

  // start optimization
  var::parm_vector vparms(num_parms_);
  var::parm_vector vparms_start(num_parms_);
  RealVector grad(num_parms_);
  RealVector grad_err(num_parms_);
  RealVector search_dir(num_parms_);
  RealMatrix sr_matrix(num_parms_,num_parms_);
  RealMatrix work_mat(num_parms_,num_parms_);
  double en, en_err;
  bool silent = true;

  if (num_opt_samples_>1) {
    vparms_start = vmc.varp_values();
  }
  else {
    vparms = vmc.varp_values();
  }

  // for all samples
  all_converged_ = true;
  for (int sample=1; sample<=num_opt_samples_; ++sample) {

    // initialize quantities for this sample
    init_sample(vmc, sample);

    // reset MK statistics
    mk_statistic_.reset();
    mk_statistic_en_.reset();

    // starting parameters (for 1 sample, already initialized)
    if (num_opt_samples_>1) vparms = vparms_start;

    // iteration count
    int iter_count = 0;
    exit_status status = exit_status::notconvgd;

   /*--------------------------------------------------------------
    * Stochastic CG iteraions 
    *--------------------------------------------------------------*/
    if (cg_maxiter_>0) {
      RealVector previous_grad(num_parms_);
      RealVector stochastic_grad(num_parms_);

      //if (print_progress_) print_progress(std::cout, iter_count, "Stochastic CG");
      //if (print_log_) print_progress(logfile_, iter_count, "Stochastic CG");

      // compute initial quantities 
      vmc.run(vparms,en,en_err,grad,grad_err,silent);
      //previous_grad = grad;
      stochastic_grad = grad;
      search_dir = -grad;

      for (int cg_iter=1; cg_iter<=cg_maxiter_; ++cg_iter) {
        iter_count++;

        if (print_progress_) print_progress(std::cout, iter_count, "Stochastic CG");
        if (print_log_) print_progress(logfile_, iter_count, "Stochastic CG");

        // max_norm (of components not hitting boundary) 
        double gnorm = grad.squaredNorm();
        double proj_norm = projected_gnorm(vparms,grad,lbound_,ubound_);
        gnorm = std::min(gnorm, proj_norm);

        // MK test: add data to Mann-Kendall statistic
        mk_statistic_ << vparms;
        mk_statistic_en_ << en;
        iter_energy_.push_back(en);
        iter_energy_err_.push_back(en_err);
        iter_gnorm_.push_back(gnorm);
        if (mk_statistic_.is_full()) {
          iter_energy_.pop_front();
          iter_energy_err_.pop_front();
          iter_gnorm_.pop_front();
        }
        // series average of gnorm
        double avg_gnorm = series_avg(iter_gnorm_);
        // MK trend
        int trend_elem;
        double mk_trend = mk_statistic_.elem_max_trend(trend_elem);
        double en_trend = mk_statistic_en_.max_trend();

        // print progress 
        if (print_progress_) {
          print_progress(std::cout,vparms,en,en_err,grad,search_dir,
            gnorm,avg_gnorm,en_trend,mk_trend,trend_elem);
        }
        if (print_log_) {
          print_progress(logfile_,vparms,en,en_err,grad,search_dir,
            gnorm,avg_gnorm,en_trend,mk_trend,trend_elem);
        }
        // file outs
        file_energy_<<std::setw(6)<<iter_count<<std::scientific<<std::setw(16)<<en; 
        file_energy_<<std::fixed<<std::setw(10)<<en_err<<std::endl<<std::flush;
        file_vparms_<<std::setw(6)<<iter_count; 
        for (int i=0; i<vparms.size(); ++i) file_vparms_<<std::setw(15)<<vparms[i];
        file_vparms_<<std::endl<<std::flush;

       /*--------------------------------------------------------------
        * Convergence criteria
        *--------------------------------------------------------------*/
        if (mk_statistic_.is_full() && en_trend<=0.1 && avg_gnorm<=grad_tol_) {
          mk_statistic_.get_series_avg(vparms);
          status = exit_status::converged;
          break;
        }
        else if (mk_statistic_.is_full() && mk_trend<=mk_thresold_ 
          && avg_gnorm<=grad_tol_) {
          // converged, add data point to store
          mk_statistic_.get_series_avg(vparms);
          status = exit_status::converged;
          break;
        }
        else if (cg_iter>=cg_maxiter_ || iter_count>= maxiter_) {
          status = exit_status::maxiter;
          break;
        }
        else if (!boost::filesystem::exists(life_fname_)) {
          status = exit_status::terminated;
          break;
        }
        else {
          // continue
          if (fixedstep_iter_ % refinement_cycle_ == 0) {
            refinement_level_++;
            search_tstep_ *= 0.5;
          }
        }

       /*----------------------------------------------------------------
        * Not converged yet. 
        * Step ahead by either Line Search OR by a fixed sized step. 
        *----------------------------------------------------------------*/
        previous_grad = grad;
        if (cg_iter <= num_probls_steps_) {
          line_search(vmc,vparms,en,en_err,grad,grad_err,search_dir);
        }
        else {
          // Fixed sized step beyond this iteration
          vparms.noalias() += search_tstep_ * search_dir;
          vparms = lbound_.cwiseMax(vparms.cwiseMin(ubound_));
          vmc.run(vparms,en,en_err,grad,grad_err,silent);
          fixedstep_iter_++;
          if (print_progress_) {
            std::cout<<" line search =  constant step\n";
            std::cout<<" step size   =  "<<search_tstep_<<"\n";
          } 
          if (print_log_) {
            logfile_<<" line search =  constant step\n";
            logfile_<<" step size   =  "<<search_tstep_<<"\n";
          }
          // run vmc at updated parameters
          vmc.run(vparms,en,en_err,grad,grad_err,silent);
        }
       /*----------------------------------------------------------------
        * New search direction from Stochastic CG
        *----------------------------------------------------------------*/
        stochastic_CG(grad, previous_grad, stochastic_grad, search_dir);
      } // CG iterations
    }

    // reset MK statistics
    mk_statistic_.reset();
    mk_statistic_en_.reset();

   /*--------------------------------------------------------------
    * Stochastic Reconfiguration iteraions 
    *--------------------------------------------------------------*/
    // counter for last data points around minimum for final average
    bool conv_condition_reached = false;
    int final_iter_count = 0;

    //exit_status status = exit_status::notconvgd;
    for (int sr_iter=1; sr_iter<=maxiter_; ++sr_iter) {
      iter_count++;

      // first message 
      if (print_progress_) print_progress(std::cout, iter_count, "SR");
      if (print_log_) print_progress(logfile_, iter_count, "SR");

      // SR search direction
      vmc.run(vparms,en,en_err,grad,grad_err,sr_matrix,silent);
      stochastic_reconf(grad,sr_matrix,work_mat,search_dir);

      // max_norm (of components not hitting boundary) 
      double gnorm = grad.squaredNorm();
      double proj_norm = projected_gnorm(vparms,grad,lbound_,ubound_);
      gnorm = std::min(gnorm, proj_norm);

      // MK test: add data to Mann-Kendall statistic
      mk_statistic_ << vparms;
      mk_statistic_en_ << en;
      iter_energy_.push_back(en);
      iter_energy_err_.push_back(en_err);
      iter_gnorm_.push_back(gnorm);
      if (mk_statistic_.is_full()) {
        iter_energy_.pop_front();
        iter_energy_err_.pop_front();
        iter_gnorm_.pop_front();
      }
      // series average of gnorm
      double avg_gnorm = series_avg(iter_gnorm_);
      // MK trend
      int trend_elem;
      double mk_trend = mk_statistic_.elem_max_trend(trend_elem);
      double en_trend = mk_statistic_en_.max_trend();

      // print progress 
      if (print_progress_) {
        print_progress(std::cout,vparms,en,en_err,grad,search_dir,
          gnorm,avg_gnorm,en_trend,mk_trend,trend_elem);
      }
      if (print_log_) {
        print_progress(logfile_,vparms,en,en_err,grad,search_dir,
          gnorm,avg_gnorm,en_trend,mk_trend,trend_elem);
      }

      // file outs
      file_energy_<<std::setw(6)<<iter_count<<std::scientific<<std::setw(16)<<en; 
      file_energy_<<std::fixed<<std::setw(10)<<en_err<<std::endl<<std::flush;
      file_vparms_<<std::setw(6)<<iter_count; 
      for (int i=0; i<vparms.size(); ++i) file_vparms_<<std::setw(15)<<vparms[i];
      file_vparms_<<std::endl<<std::flush;


     /*--------------------------------------------------------------
      * Convergence criteria
      *--------------------------------------------------------------*/
      if (conv_condition_reached) {
        if (print_progress_) {
          std::cout<<" final iter  =  "<<final_iter_count<<"/"<<mk_series_len_<<"\n";
        } 
        if (print_log_) {
          logfile_<<" final iter  =  "<<final_iter_count<<"/"<<mk_series_len_<<"\n";
        }
        final_iter_count++;
        if (final_iter_count >= mk_series_len_) {
          mk_statistic_.get_series_avg(vparms);
          status = exit_status::converged;
          break;
        }
      }

      if (mk_statistic_.is_full() && en_trend<=0.1 && avg_gnorm<=grad_tol_) {
        conv_condition_reached = true;
        //mk_statistic_.get_series_avg(vparms);
        //status = exit_status::converged;
        //break;
      }
      else if (mk_statistic_.is_full() && mk_trend<=mk_thresold_ 
        && avg_gnorm<=grad_tol_) {
        conv_condition_reached = true;
        // converged, add data point to store
        //mk_statistic_.get_series_avg(vparms);
        //status = exit_status::converged;
        //break;
      }
      else if (sr_iter>=maxiter_ || iter_count>= maxiter_) {
        status = exit_status::maxiter;
        break;
      }
      else if (!boost::filesystem::exists(life_fname_)) {
        status = exit_status::terminated;
        break;
      }
      else {
        // continue
        if (fixedstep_iter_ % refinement_cycle_ == 0) {
          refinement_level_++;
          search_tstep_ *= 0.5;
        }
      }

      /*----------------------------------------------------------------
       * Not converged yet. 
       * Step ahead by either Line Search OR by a fixed sized step. 
       *----------------------------------------------------------------*/
      if (sr_iter <= num_probls_steps_) {
        // Line search
        line_search(vmc,vparms,en,en_err,grad,grad_err,search_dir);
      }
      else {
        // Fixed sized step beyond this iteration
        vparms.noalias() += search_tstep_ * search_dir;
        vparms = lbound_.cwiseMax(vparms.cwiseMin(ubound_));
        fixedstep_iter_++;
        if (print_progress_) {
          std::cout<<" line search =  constant step\n";
          std::cout<<" step size   =  "<<search_tstep_<<"\n";
        } 
        if (print_log_) {
          logfile_<<" line search =  constant step\n";
          logfile_<<" step size   =  "<<search_tstep_<<"\n";
        }
      }
    } // iterations

    // Iterations over for this sample, finalize
    if (print_progress_) {
      std::cout << "#" << std::string(60, '-') << std::endl;
      if (status == exit_status::converged) {
        std::cout<<" Converged!"<<std::endl;
      }
      else if (status == exit_status::maxiter) {
        std::cout<<" Iterations exceeded (NOT converged)"<<std::endl;
      }
      else if (status == exit_status::terminated) {
        std::cout<<" Iterations terminated (NOT converged)"<<std::endl;
      }
      else {
        std::cout <<" NOT converged"<<std::endl;
      }
      std::cout << "#" << std::string(60, '-') << std::endl << std::flush;
    }
    if (print_log_) {
      logfile_ << "#" << std::string(60, '-') << std::endl;
      if (status == exit_status::converged) {
        logfile_<<" Converged!"<<std::endl;
      }
      else if (status == exit_status::maxiter) {
        logfile_<<" Iterations exceeded (NOT converged)"<<std::endl;
      }
      else if (status == exit_status::terminated) {
        logfile_<<" Iterations terminated (NOT converged)"<<std::endl;
      }
      else {
        logfile_ <<" NOT converged"<<std::endl;
      }
      logfile_ << "#" << std::string(60, '-')<<std::endl<<std::endl<< std::flush;
    }
    if (status != exit_status::converged) {
      all_converged_ = false;
    }
    finalize_sample(status);

  }// end of samples

  return 0;
}


int Optimizer::finalize_sample(const exit_status& status)
{
  file_energy_.close();
  file_vparms_.close();
  for (int i=0; i<iter_energy_.size(); ++i) {
    optimal_energy_ << iter_energy_[i];
    energy_error_bar_ << iter_energy_err_[i];
  }
  for (const auto& p : mk_statistic_.data_series()) {
    optimal_parms_ << p;
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
  optimal_energy_.fs()<<std::setw(10)<<num_vmc_samples_; 
  if (status == exit_status::converged) {
    optimal_energy_.fs()<<std::setw(11)<<"CONVERGED"<<std::setw(7)<<0<<std::endl;
  }
  else {
    optimal_energy_.fs()<<std::setw(11)<<"NOT_CONVD"<<std::setw(7)<<0<<std::endl;
  }
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
  optimal_energy_.fs()<<std::setw(10)<<num_vmc_samples_; 
  if (all_converged_) 
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

  // reset for the next sample
  optimal_energy_.reset();
  energy_error_bar_.reset();
  optimal_parms_.reset();
  return 0;
}

std::ostream& Optimizer::print_progress(std::ostream& os, 
  const int& iter, const std::string& algorithm)
{
  std::ios state(NULL);
  state.copyfmt(os);
  os << "#" << std::string(60, '-') << std::endl;
  os << std::left;
  os << " iteration   =  "<<std::setw(6)<<iter<<"\n";
  os << " step        =  "<< algorithm <<"\n";
  os.copyfmt(state);
  return os;
}

std::ostream& Optimizer::print_progress(std::ostream& os, 
  const var::parm_vector& vparms, const double& energy,
  const double& error_bar, const RealVector& grad, const RealVector& search_dir,
  const double& gnorm, const double& avg_gnorm, const double& en_trend, 
  const double& mk_trend, const int& trend_elem) 
{
  std::ios state(NULL);
  state.copyfmt(os);
  os <<std::scientific<<std::uppercase<<std::setprecision(6)<<std::right;
  os << " energy      ="<<std::setw(14)<<energy << "   (+/-) ";
  os <<std::fixed<<std::setw(10)<<error_bar<<"\n";
  os <<std::scientific<<std::right;
  os <<std::scientific<<std::uppercase<<std::setprecision(4)<<std::right;
  os << " varp        =";
  for (int i=0; i<num_parms_print_; ++i) os<<std::setw(12)<<vparms[i];
  os << "\n";
  os << " grad        =";
  for (int i=0; i<num_parms_print_; ++i) os<<std::setw(12)<<grad[i];
  os << "\n";
  os << " search_dir  =";
  for (int i=0; i<num_parms_print_; ++i) os<<std::setw(12)<<search_dir[i];
  os << "\n";
  os << " gnorm       ="<<std::setw(12)<< gnorm;
  os << " (avg ="<<std::setw(11)<<avg_gnorm<<")\n";
  os <<std::fixed<<std::setprecision(6)<<std::right;
  os <<" MK trend    = "<<std::setw(9)<<en_trend<<"  ";
  os <<std::setw(9)<<mk_trend;
  os << " ("<<trend_elem<<")"<<"\n"; 
  os.copyfmt(state);
  return os;
}


/*
int Optimizer::genetic_algorthim(const RealVector& grad, RealMatrix& srmat, 
  RealMatrix& srmat_inv, RealVector& search_dir)
{
  while (not converged) {
    vmc.run(vparms,en,en_err,true);
  }
}
*/


int Optimizer::stochastic_reconf(const RealVector& grad, RealMatrix& srmat, 
  RealMatrix& srmat_inv, RealVector& search_dir)
{

  // conditioning
  for (int i=0; i<num_parms_; ++i) srmat(i,i) *= (1.0+stabilizer_);

  // search dir after cutting off redundant directions
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(srmat, Eigen::ComputeThinU | Eigen::ComputeThinV); 
  double lmax = svd.singularValues()(0);
  for (int i=0; i<num_parms_; ++i) {
    for (int j=0; j<num_parms_; ++j) {
      double sum = 0.0;
      for (int m=0; m<num_parms_; ++m) {
        double sval = svd.singularValues()(m);
        if (sval/lmax < w_svd_cut_) break;
        sum += svd.matrixV()(i,m)*svd.matrixU()(j,m)/sval;
      }
      srmat_inv(i,j) = sum;
    }
  }
  search_dir = -srmat_inv*grad;

  return 0;
}

double Optimizer::line_search(VMCRun& vmc, RealVector& vparms, 
  double& en, double& en_err, RealVector& grad, RealVector& grad_err, 
  const RealVector& search_dir)
{
  if (print_progress_) std::cout << " line search =  ";
  if (print_log_) logfile_ << " line search =  "; 
  std::string message;
  bool search_done = false;
  int search_step = 0;
  double net_dt = 0.0;
  // start LS
  double dt = probLS_.start(en,en_err,grad,grad_err,search_dir);
  net_dt += dt;

  // do steps
  while (true) {
    search_step++;
    if (print_progress_) {
      std::cout << search_step << ".." << std::flush;
    }
    if (print_log_) {
      logfile_ << search_step << ".."; //<< std::flush; 
    }

    // run with next parameters
    vparms += dt*search_dir;
    vparms = lbound_.cwiseMax(vparms.cwiseMin(ubound_));
    vmc.run(vparms,en,en_err,grad,grad_err,true);

    // find next step size by probabilistic line search method
    dt = probLS_.do_step(en,en_err,grad,grad_err,search_dir,search_done,message);
    net_dt += dt;

    // check if done
    if (search_done) {
      if (print_progress_) {
        //std::cout<<" done\n";
        //std::cout<<" step size   =  "<<net_dt<<" ("<<message<<")\n";
        std::cout<<" "<<message<<"\n";
        std::cout<<" step size   =  "<<net_dt<<"\n";
      } 
      if (print_log_) {
        logfile_<<" done\n";
        logfile_<<" step size   =  "<<net_dt<<" ("<<message<<")\n";
      }
      break;
    } 
    //std::cout << "dt = " << dt << "\n";
  }
  // return step size
  return net_dt;
}


double Optimizer::projected_gnorm(const RealVector& x, const RealVector& grad, 
  const RealVector& lb, const RealVector& ub) const
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

double Optimizer::series_avg(const std::deque<double>& data) const
{
  if (data.size() == 0) return 0.0;
  double sum = 0.0;
  for (const auto& d : data) sum += d;
  return sum/data.size();
}

int Optimizer::stochastic_CG(const RealVector& grad, const RealVector& grad_xprev, 
  RealVector& stochastic_grad, RealVector& search_dir)
{
  // Norm of previous stochastic gradient
  double gnorm_prev = stochastic_grad.squaredNorm();

  // Update stocahastic gradient (SARAH algorithm)
  RealVector grad_diff = grad - grad_xprev;
  stochastic_grad += grad_diff;

  // Conjugate search direction
  // compute beta
  double beta;
  if (cg_algorithm_==cg_type::FR) {
    beta = stochastic_grad.squaredNorm()/gnorm_prev;
  }
  else if (cg_algorithm_==cg_type::PR) {
    beta = stochastic_grad.dot(grad_diff)/gnorm_prev;
  }
  else if (cg_algorithm_==cg_type::HS) {
    beta = stochastic_grad.dot(grad_diff)/search_dir.dot(grad_diff);
  }
  else if (cg_algorithm_==cg_type::DY) {
    beta = stochastic_grad.squaredNorm()/search_dir.dot(grad_diff);
  }
  else {
    throw std::range_error("Optimizer::stochastic_CG: unknown CG method");
  }
  //std::cout << "beta = " << beta << "\n"; //getchar();

  // new search direction
  search_dir = -stochastic_grad+beta*search_dir;

  return 0;
}


} // end namespace vmc














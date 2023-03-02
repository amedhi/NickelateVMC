/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-12 13:20:56
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-17 07:24:49
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <boost/filesystem.hpp>
#include <Eigen/SVD>
#include "vmc.h"

namespace vmc {

VMC::VMC(const input::Parameters& inputs) 
  : lattice(inputs) 
  , model(inputs,lattice) 
  , config(inputs,lattice,model)
  , num_varparms_{config.num_varparms()}
{
  // seed random generator
  rng_seed_ = inputs.set_value("rng_seed", 1);
  config.rng().seed(rng_seed_);
  // mc parameters
  num_sites_ = lattice.num_sites();
  num_measure_steps_ = inputs.set_value("measure_steps", 0);
  num_warmup_steps_ = inputs.set_value("warmup_steps", 0);
  min_interval_ = inputs.set_value("min_interval", 0);
  max_interval_ = inputs.set_value("max_interval", 0);
  check_interval_ = std::max(1,num_measure_steps_/10);

  // disorder
  if (site_disorder_.check(inputs)) {
    site_disorder_.init(inputs,lattice,model,config,config.rng());
    model.add_disorder_term("disorder", model::op::ni_sigma());
    model.finalize(lattice);
  }
  //if (model.have_disorder_term()) 

  // output directory
  std::string outdir = inputs.set_value("prefix", "");
  prefix_ = "./"+outdir+"/";
  boost::filesystem::path prefix_dir(prefix_);
  boost::filesystem::create_directory(prefix_dir);

  // x-variables (assuming max 1)
  xvar_names_.resize(1);
  xvar_values_.resize(1);
  xvar_names_[0]=inputs.set_value("as_func_of","hole_doping");
  if (xvar_names_[0]=="hole_doping")
    xvar_values_[0] = config.hole_doping();
  else {
    xvar_values_[0] = inputs.set_value(xvar_names_[0],0.0);
  }

  // observables
  make_info_str(inputs);
  observables.headstream() << info_str_.str();
  observables.init(inputs,lattice,model,config,prefix_,num_measure_steps_);
  observables.as_functions_of(xvar_names_);
}

void VMC::seed_rng(const int& seed)
{
  config.rng().seed(seed);
}

int VMC::disorder_start(const input::Parameters& inputs, 
  const unsigned& disorder_config, const run_mode& mode, const bool& silent)
{
  site_disorder_.set_current_config(disorder_config);
  //site_disorder_.save_optimal_parms(config.vparm_values());
  run_mode_ = mode;
  bool with_psi_grad = false;
  switch (mode) {
    case run_mode::normal:
      if (observables.energy_grad()) with_psi_grad = true;
      break;
    case run_mode::en_function:
      observables.switch_off();
      observables.energy().setup(lattice,model);
      observables.energy_grad().setup(config);
      break;
    case run_mode::sr_function:
      observables.switch_off();
      observables.energy().setup(lattice,model);
      observables.energy_grad().setup(config);
      observables.sr_matrix().setup(lattice,config);
      break;
  }
  silent_mode_ = silent;
  if (inputs.have_option_quiet()) silent_mode_ = true;
  //site_disorder_.get_optimal_parms();
  //return config.build(lattice, inputs, with_psi_grad);
  return config.build(lattice,site_disorder_.get_optimal_parms(),with_psi_grad);
}

int VMC::set_run_mode(const run_mode& mode)
{
  run_mode_ = mode;
  switch (mode) {
    case run_mode::normal:
      break;
    case run_mode::en_function:
      observables.switch_off();
      observables.energy().setup(lattice,model,true); // total energy
      observables.energy_grad().setup(config);
      break;
    case run_mode::sr_function:
      observables.switch_off();
      observables.energy().setup(lattice,model,true); // total energy
      observables.energy_grad().setup(config);
      observables.sr_matrix().setup(lattice,config);
      break;
  }
  return 0;
}

int VMC::start(const input::Parameters& inputs, const run_mode& mode, 
  const bool& silent)
{
  config.rng().seed(rng_seed_);
  silent_mode_ = silent;
  if (inputs.have_option_quiet()) silent_mode_ = true;
  set_run_mode(mode);

  // build config
  bool with_psi_grad = false;
  if (observables.energy_grad()) with_psi_grad = true;
  config.build(lattice, inputs, with_psi_grad);

  // xvariables
  if (xvar_names_[0] == "hole_doping") {
    xvar_values_[0] = config.hole_doping();
  }
  else {
    xvar_values_[0] = inputs.set_value(xvar_names_[0],0.0);
  }
  return 0;
}

int VMC::start(const var::parm_vector& varp, const run_mode& mode, const bool& silent) 
{
  silent_mode_ = silent;
  set_run_mode(mode);
  config.rng().seed(rng_seed_);
  bool with_psi_grad = false;
  if (observables.energy_grad()) with_psi_grad = true;
  return config.build(lattice,varp,with_psi_grad);
  return 0;
}

int VMC::build_config(const Eigen::VectorXd& varp, const bool& with_psi_grad) 
{
  config.rng().seed(rng_seed_);
  return config.build(lattice, varp, with_psi_grad);
}

int VMC::sr_function(const var::parm_vector& varp, double& en_mean, 
  double& en_err, RealVector& grad, RealVector& grad_err, RealMatrix& sr_matrix,
  const bool& silent)
{
  config.rng().seed(rng_seed_);
  silent_mode_ = silent;
  set_run_mode(run_mode::sr_function);
  // build the config from the variational parameters
  bool with_psi_grad = true;
  config.build(lattice, varp, with_psi_grad);
  // run the simulation
  run_simulation();
  // energy
  en_mean = observables.energy().mean(0);
  en_err = observables.energy().stddev(0);
  // gradient
  grad = observables.energy_grad().mean_data();
  grad_err = observables.energy_grad().stddev_data();
  // sr matrix
  observables.sr_matrix().get_matrix(sr_matrix);
  return 0;
}

int VMC::en_function(const var::parm_vector& varp, double& en_mean, 
  double& en_stddev, RealVector& grad, RealVector& grad_stddev,
  const bool& silent)
{
  config.rng().seed(rng_seed_);
  silent_mode_ = silent;
  set_run_mode(run_mode::en_function);
  // build the config from the variational parameters
  bool with_psi_grad = true;
  config.build(lattice, varp, with_psi_grad);
  // run the simulation
  run_simulation();
  // results
  en_mean = observables.energy().mean(0);
  en_stddev = observables.energy().stddev(0);
  // gradient
  grad = observables.energy_grad().mean_data();
  grad_stddev = observables.energy_grad().stddev_data();
  return 0;
}

// simulation after optimization
/*
int VMC::run_simulation(const Eigen::VectorXd& varp)
{
  observables.switch_off();
  observables.energy().switch_on();
  config.build(lattice, varp);
  run_simulation();
  return 0;
}
*/

int VMC::run_simulation(void)
{
  std::vector<int> empty_bc_list;
  return run_simulation(num_measure_steps_,empty_bc_list);
}

int VMC::run_simulation(const int& sample_size)
{
  std::vector<int> empty_bc_list;
  return run_simulation(sample_size,empty_bc_list);
}

int VMC::run_simulation(const std::vector<int>& bc_list)
{
  return run_simulation(num_measure_steps_,bc_list);
}

int VMC::run_simulation(const int& sample_size, const std::vector<int>& bc_list)
{
  /*--------------------------------------------------------------
   * Assumming everything ready, run actual simulation here
   *--------------------------------------------------------------*/
  assert(sample_size>=0);

  // update 'sample_size' dependent things
  check_interval_ = std::max(1,sample_size/10);
  observables.reset_batch_limit(sample_size);

  // Take care of BC twists
  std::vector<int> bc_twists;
  if (bc_list.size()==0) {
    for (int bc=0; bc<lattice.num_boundary_twists(); ++bc) {
      bc_twists.push_back(bc);
    }
  }
  else {
    bc_twists = bc_list;
    //for (const int& bc : bc_list) bc_twists.push_back(bc);
  }

  /*******************************************
   * Only one and default BC
   * *****************************************/
  if (bc_twists.size()==1 && bc_twists[0]==0) {
    // initialize
    observables.reset();

    // warmup run
    if (!silent_mode_) std::cout << " warming up... " << std::flush;
    do_warmup_run();
    if (!silent_mode_) std::cout << "done\n" << std::flush;
    // measuring run
    do_measure_run(sample_size);
    // finalize
    observables.finalize();
    if (!silent_mode_) {
      std::cout << " simulation done\n";
      config.print_stats();
    }
  }
  /*******************************************
   * Average over several BCs
   * *****************************************/
  else {
    // initialize
    observables.reset();

    for (const auto& bc: bc_twists) {
      if (!silent_mode_) {
        std::cout << "\n-------------------------------------" << std::endl;
        std::cout << " Running for BC twist - "<<bc+1 << " / ";
        std::cout << bc_twists.size() << std::endl;
        std::cout << "-------------------------------------\n" << std::flush;
      }
      lattice.reset_boundary_twist(bc);
      config.rng().seed(rng_seed_);
      config.rebuild(lattice);

      // warmup run
      if (!silent_mode_) std::cout << " warming up... " << std::flush;
      do_warmup_run();
      if (!silent_mode_) std::cout << "done\n" << std::flush;
      // measuring run
      do_measure_run(sample_size);
      // finalization needed only for "Energy Gradient"
      observables.finalize();

      //std::cout<<"Energy (bc="<<bc<<") = "<<observables.energy().result_str()<<"\n";
      // getchar();
    }
    //std::cout << observables.energy().result_str(-1) << "\n";
    //getchar();
    if (!silent_mode_) {
      std::cout << " simulation done\n";
      config.print_stats();
    }
  }
  return 0;
}

int VMC::reset_observables(void) 
{
  observables.reset();
  return 0;
}

int VMC::do_warmup_run(void) 
{
  for (int n=0; n<num_warmup_steps_; ++n) {
    config.update_state();
  }
  return 0;
}

int VMC::do_measure_run(const int& num_samples) 
{
  int skip_count = min_interval_;
  int measurement_count = 0;
  while (measurement_count < num_samples) {
    config.update_state();
    if (skip_count >= min_interval_) {
      if (config.accept_ratio()>0.5 || skip_count==max_interval_) {
        skip_count = 0;
        config.reset_accept_ratio();
        observables.do_measurement(lattice,model,config,site_disorder_);
        ++measurement_count;
        if (!silent_mode_) print_progress(measurement_count, num_samples);
      }
    }
    skip_count++;
  }
  return 0;
}

void VMC::print_results(void) 
{
  //if (observables.energy_grad()) //finalize_energy_grad();
  //observables.print(config.wavefunc().hole_doping());
  //observables.print_results(config.vparm_vector());
  //observables.print_results(config.hole_doping());
  observables.print_results(xvar_values_);
  //std::cout << observables.energy().with_statistic() << "\n";
}

void VMC::print_progress(const int& num_measurement, const int& num_measure_steps) const
{
  if (num_measurement%check_interval_==0) {
    std::ios state(NULL);
    state.copyfmt(std::cout);
    std::cout<<std::fixed<<std::setprecision(1)<<std::right;
    double w = double(100.0*num_measurement)/num_measure_steps;
    std::cout<<" measurement = "<< w <<" %\n";
    std::cout.copyfmt(state);
  }
}

void VMC::make_info_str(const input::Parameters& inputs)
{
  info_str_.clear();
  copyright_msg(info_str_);
  info_str_ << "# "<< inputs.job_id() <<"\n"; 
  info_str_ << model.info_str(); 
  info_str_ << config.info_str(); 
  info_str_ << "# Samples = " << num_measure_steps_;
  info_str_ << ", warmup = " << num_warmup_steps_;
  info_str_ << ", min_interval = " << min_interval_;
  info_str_ << ", max_interval = " << max_interval_ << "\n";
}

void VMC::copyright_msg(std::ostream& os)
{
  os << "#" << std::string(72,'-') << "\n";
  os << "#" << " Program: VMC Simulation\n";
  os << "#" << "          (c) Amal Medhi <amedhi@iisertvm.ac.in>\n";
  os << "#" << std::string(72,'-') << "\n";
}

/*
* The OLD one
int VMC::run_simulation(const int& sample_size, const std::vector<int>& bc_list)
{
  // take care of default argument
  std::vector<int> bc_twists;
  if (bc_list.size()==1 && bc_list[0]==-1) {
    for (int bc=0; bc<lattice.num_boundary_twists(); ++bc) {
      bc_twists.push_back(bc);
    }
  }
  else {
    bc_twists = bc_list;
    //for (const int& bc : bc_list) bc_twists.push_back(bc);
  }

  // *******************************************
  // * Only one and default BC
  // *******************************************
  if (bc_twists.size()==1 && bc_twists[0]==0) {
    // initialize
    observables.reset();
    int num_measure_steps = num_measure_steps_;
    if (sample_size>0) num_measure_steps = sample_size;

    // warmup run
    if (!silent_mode_) std::cout << " warming up... " << std::flush;
    for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
    if (!silent_mode_) std::cout << "done\n" << std::flush;
    // measuring run
    int skip_count = min_interval_;
    int measurement_count = 0;
    while (measurement_count < num_measure_steps) {
      config.update_state();
      if (skip_count >= min_interval_) {
        if (config.accept_ratio()>0.5 || skip_count==max_interval_) {
          skip_count = 0;
          config.reset_accept_ratio();
          observables.do_measurement(lattice,model,config,site_disorder_);
          ++measurement_count;
          if (!silent_mode_) print_progress(measurement_count, num_measure_steps);
        }
      }
      skip_count++;
    }
    // finalize
    observables.finalize();
    if (!silent_mode_) {
      std::cout << " simulation done\n";
      config.print_stats();
    }
  }
  // *******************************************
  // Average over several BCs
  // *******************************************
  else {
    // initialize
    int num_measure_steps = num_measure_steps_;
    if (sample_size>0) num_measure_steps = sample_size;
    observables.reset();
    for (const auto& bc: bc_twists) {
      if (!silent_mode_) {
        std::cout << "\n-------------------------------------" << std::endl;
        std::cout << " Running for BC twist - "<<bc+1 << " / ";
        std::cout << bc_twists.size() << std::endl;
        std::cout << "-------------------------------------\n" << std::flush;
      }
      lattice.reset_boundary_twist(bc);
      config.rng().seed(rng_seed_);
      config.rebuild(reset;
      //observables.reset();

      // warmup run
      if (!silent_mode_) std::cout << " warming up... " << std::flush;
      for (int n=0; n<num_warmup_steps_; ++n) config.update_state();
      if (!silent_mode_) std::cout << "done\n" << std::flush;
      // measuring run
      int skip_count = min_interval_;
      int measurement_count = 0;
      while (measurement_count < num_measure_steps) {
        config.update_state();
        if (skip_count >= min_interval_) {
          if (config.accept_ratio()>0.5 || skip_count==max_interval_) {
            skip_count = 0;
            config.reset_accept_ratio();
            observables.do_measurement(lattice,model,config,site_disorder_);
            ++measurement_count;
            if (!silent_mode_) print_progress(measurement_count, num_measure_steps);
          }
        }
        skip_count++;
      }
      // finalization needed only for "Energy Gradient"
      observables.finalize();
      std::cout << observables.energy().result_str(-1) << "\n";
    }
    //std::cout << observables.energy().result_str(-1) << "\n";
    //getchar();
    if (!silent_mode_) {
      std::cout << " simulation done\n";
      config.print_stats();
    }
  }
  return 0;
}
*/



} // end namespace vmc






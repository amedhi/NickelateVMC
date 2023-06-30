/*---------------------------------------------------------------------------
* @Author: Amal Medhi
* @Date:   2022-10-15 14:57:27
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2023-06-26 13:19:30
* Copyright (C) 2015-2022 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
*----------------------------------------------------------------------------*/
#include "./vmcrun.h"

namespace vmc {

VMCRun::VMCRun(const input::Parameters& inputs) 
  : VMC(inputs)
{
}

// Serial run, normal mode
int VMCRun::run(const input::Parameters& inputs)
{
  if (!inputs.have_option_quiet()) std::cout << " starting vmc run\n";
  start(inputs, run_mode::normal);
  run_simulation();
  print_results();
  return 0;
}

// MPI run, normal mode
int VMCRun::run(const input::Parameters& inputs, const mpi::mpi_communicator& mpi_comm) 
{
  if (mpi_comm.is_master() && !inputs.have_option_quiet()) {
  	if (mpi_mode_ == mpi_mode::NORMAL) {
      std::cout << " starting vmc run in MPI_mode::NORMAL\n";
  	}
  	if (mpi_mode_ == mpi_mode::BC_TWIST) {
      std::cout << " starting vmc run in MPI_mode::BC_TWIST\n";
  	}
  }

  bool silent;
  if (mpi_comm.is_master() && !inputs.have_option_quiet()) {
    std::cout << " *progress in master*\n";
    silent = false;
    start(inputs,run_mode::normal,silent);
  }
  else {
    silent = true;
    start(inputs,run_mode::normal,silent);
  }

  if (mpi_mode_ == mpi_mode::NORMAL) {
    //std::cout << " running "<<measure_samples_<<" samples in p"
    //	<<mpi_comm.rank()<<"\n";
  	run_simulation(measure_samples_);
  }
  if (mpi_mode_ == mpi_mode::BC_TWIST) {
    //std::cout << " running BC_TWISTS: "<<bc_list_.front()+1<<"-"
    //	<<bc_list_.back()+1<<" in p"<<mpi_comm.rank()<<"\n";
  	run_simulation(bc_list_);
  }

  // collect samples
  if (mpi_comm.is_master()) {
    if (!inputs.have_option_quiet()) {
      std::cout << " collecting data from other procs...";
      std::cout << std::flush;
    }

    // Collect results: Sequential order
    /*
    for (const mpi::proc& p : worker_procs_) {
      //std::cout << "reciving results from p = " << p << "\n";
      MPI_recv_results(mpi_comm, p, mpi::MP_data_samples);
    }*/

    // Collect results: First-come-first-serve basis
    std::set<mpi::proc> working_procs; 
    for (const mpi::proc& p : worker_procs_) working_procs.insert(p);
    //std::cout << "\n";
    while (working_procs.size() > 0) {
      mpi::mpi_status msg = mpi_comm.probe(mpi::any_source, mpi::MP_data_samples);
      //std::cout << "reciving results from p = " << msg.source() << "\n";
      MPI_recv_results(mpi_comm, msg.source(), mpi::MP_data_samples);
      working_procs.erase(msg.source());
    }

  }
  else {
    //std::cout << "sending results from p = " << mpi_comm.rank() << "\n";
    MPI_send_results(mpi_comm, mpi_comm.master(), mpi::MP_data_samples);
  }

  // finalize results
  if (mpi_comm.is_master()) {
    finalize_results();
  	print_results();
    if (!inputs.have_option_quiet()) {
      std::cout << " done\n";
      std::cout << std::flush;
    }
  }

  // for safetly, all together here
  mpi_comm.barrier(); 

  return 0;
}

// Dynamic run, sr_function mode
int VMCRun::run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err, RealMatrix& sr_matrix,
	const bool& silent_mode)
{
  if (!mpi_comm_.is_master()) {
    throw std::logic_error("VMCRun::run: this should be the master");
  }
  // if MPI
  if (mpi_comm_.size()>1) {
  	return master_run(vparms,en,en_err,grad,grad_err,sr_matrix,silent_mode);
  }
  else {
  	// if serial
  	return sr_function(vparms,en,en_err,grad,grad_err,sr_matrix,silent_mode);
  }
}

// master controlling, sr_function mode run
int VMCRun::master_run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err, RealMatrix& sr_matrix, 
	const bool& silent_mode)
{
  // start worker simulation
  for (const mpi::proc& p : worker_procs_) {
    mpi_comm_.send(p, mpi::MP_variational_parms, vparms);
    mpi_comm_.send(p, mpi::MP_start_simulation, run_mode::sr_function);
  }
  start(vparms,run_mode::sr_function,silent_mode);
  if (mpi_mode_ == mpi_mode::NORMAL) {
  	run_simulation(measure_samples_);
  }
  if (mpi_mode_ == mpi_mode::BC_TWIST) {
  	run_simulation(bc_list_);
  }

  // Collect results: Sequential order
  /*
  for (const mpi::proc& p : worker_procs_) {
    MPI_recv_results(mpi_comm_,p,mpi::MP_data_samples);
  }*/

  // Collect results: First-come-first-serve basis
  std::set<mpi::proc> working_procs; 
  for (const mpi::proc& p : worker_procs_) working_procs.insert(p);
  while (working_procs.size() > 0) {
    mpi::mpi_status msg = mpi_comm_.probe(mpi::any_source, mpi::MP_data_samples);
    MPI_recv_results(mpi_comm_, msg.source(), mpi::MP_data_samples);
    working_procs.erase(msg.source());
  }

  finalize_results();
  //print_results();

  // computed quantities
  en = observable().energy().mean(0);
  //std::cout << "energy = " << en <<"\n"; getchar();
  en_err = observable().energy().stddev(0);
  grad = observable().energy_grad().mean_data();
  grad_err = observable().energy_grad().stddev_data();
  observable().sr_matrix().get_matrix(sr_matrix);
  //std::cout << grad_err.transpose()<<"\n";

  return 0;
}

// Dynamic run, en_function mode
int VMCRun::run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err,
	const bool& silent_mode)
{
  if (!mpi_comm_.is_master()) {
    throw std::logic_error("VMCRun::run: this should be the master");
  }
  // if MPI
  if (mpi_comm_.size()>1) {
  	return master_run(vparms,en,en_err,grad,grad_err,silent_mode);
  }
  else {
  	// if serial
  	return en_function(vparms,en,en_err,grad,grad_err,silent_mode);
  }
}

// master controlling, en_function mode run
int VMCRun::master_run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err,const bool& silent_mode)
{
  // start worker simulation
  for (const mpi::proc& p : worker_procs_) {
    mpi_comm_.send(p, mpi::MP_variational_parms, vparms);
    mpi_comm_.send(p, mpi::MP_start_simulation, run_mode::en_function);
  }
  start(vparms,run_mode::en_function,silent_mode);
  if (mpi_mode_ == mpi_mode::NORMAL) {
  	run_simulation(measure_samples_);
  }
  if (mpi_mode_ == mpi_mode::BC_TWIST) {
  	run_simulation(bc_list_);
  }

  // Collect results: Sequential order
  /*for (const mpi::proc& p : worker_procs_) {
    MPI_recv_results(mpi_comm_,p,mpi::MP_data_samples);
  }*/

  // Collect results: First-come-first-serve basis
  std::set<mpi::proc> working_procs; 
  for (const mpi::proc& p : worker_procs_) working_procs.insert(p);
  while (working_procs.size() > 0) {
    mpi::mpi_status msg = mpi_comm_.probe(mpi::any_source, mpi::MP_data_samples);
    MPI_recv_results(mpi_comm_, msg.source(), mpi::MP_data_samples);
    working_procs.erase(msg.source());
  }
  // finalize results
  finalize_results();
  print_results();

  // computed quantities
  en = observable().energy().mean(0);
  en_err = observable().energy().stddev(0);
  grad = observable().energy_grad().mean_data();
  grad_err = observable().energy_grad().stddev_data();
  //std::cout << "grad_err = " << grad_err.transpose()<<"\n";
  return 0;
}


int VMCRun::start_worker_run(void)
{
  run_mode runmode=run_mode::normal;
  bool silent_mode = true;
  var::parm_vector vparms(num_varp());
  bool not_done = true;
  while (not_done) {
  	mpi::mpi_status msg = mpi_comm_.probe();
  	switch (msg.tag()) {
  	  case mpi::MP_variational_parms:
        mpi_comm_.recv(msg.source(),msg.tag(),vparms);
        break;
  	  case mpi::MP_start_simulation:
        mpi_comm_.recv(msg.source(),msg.tag(),runmode);
  		  start(vparms,runmode,silent_mode);
  		  if (mpi_mode_ == mpi_mode::NORMAL) {
  	      run_simulation(measure_samples_);
        }
        if (mpi_mode_ == mpi_mode::BC_TWIST) {
  	      run_simulation(bc_list_);
        }
    	 MPI_send_results(mpi_comm_, mpi_comm_.master(), mpi::MP_data_samples);
  	  	break;
  	  case mpi::MP_stop_simulation:
        mpi_comm_.recv(msg.source(), msg.tag());
        not_done = false;
        break;
  	  default: 
        std::cout << " >> VMCRun::start_worker_run: unexpected message\n";
  	  	break;
  	}
  }
  return 0;
}

int VMCRun::stop_worker_run(void)
{
  for (const mpi::proc& p : worker_procs_) {
    mpi_comm_.send(p, mpi::MP_stop_simulation);
  }
  return 0;
}

// for external optimizers
/*
bool VMCRun::Evaluate(const double* parameters, double* cost, double* gradient)
{
  double x = parameters[0];
  double y = parameters[1];
  cost[0] = (1.0 - x) * (1.0 - x) + 100.0 * (y - x * x) * (y - x * x);
  if (gradient) {
     gradient[0] = -2.0 * (1.0 - x) - 200.0 * (y - x * x) * 2.0 * x;
     gradient[1] = 200.0 * (y - x * x);
  }
  return true;
}

int VMCRun::NumParameters(void) const 
{ 
  return 2; 
}
*/



} // end namespace vmc
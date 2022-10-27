/*---------------------------------------------------------------------------
* @Author: Amal Medhi
* @Date:   2022-10-15 14:57:27
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2022-10-26 15:53:29
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

  // every procs
  start(inputs,run_mode::normal,true);
  if (mpi_mode_ == mpi_mode::NORMAL) {
    //std::cout << " running "<<measure_samples_<<" samples in p"
    //	<<mpi_comm.rank()<<"\n";
  	run_simulation(measure_samples_);
  }
  if (mpi_mode_ == mpi_mode::BC_TWIST) {
    //std::cout << " running BC_TWISTS: "<<bc_list_.front()+1<<"-"
    //	<<bc_list_.back()+1<<" in p"<<mpi_comm.rank()<<"\n";
  	run_simulation(-1, bc_list_);
  }

  // collect samples
  if (mpi_comm.is_master()) {
    for (const mpi::proc& p : worker_procs_) {
      //std::cout << "reciving results from p = " << p << "\n";
      MPI_recv_results(mpi_comm, p, mpi::MP_data_samples);
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
  	run_simulation(-1, bc_list_);
  }
  for (const mpi::proc& p : worker_procs_) {
    MPI_recv_results(mpi_comm_,p,mpi::MP_data_samples);
  }
  finalize_results();
  //print_results();

  // computed quantities
  en = observable().energy().mean();
  en_err = observable().energy().stddev();
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
  	run_simulation(-1, bc_list_);
  }
  for (const mpi::proc& p : worker_procs_) {
    MPI_recv_results(mpi_comm_,p,mpi::MP_data_samples);
  }
  finalize_results();
  print_results();

  // computed quantities
  en = observable().energy().mean();
  en_err = observable().energy().stddev();
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
  	      run_simulation(-1, bc_list_);
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


#ifdef ON
int VMCRun::mpirun(const input::Parameters& inputs, const mpi::mpi_communicator& mpi_comm) 
{
  if (mpi_comm.is_master()) {
    if (!inputs.have_option_quiet()) std::cout << " starting vmc run in normal MPI mode\n";
  }
  // every procs
  start(inputs,run_mode::normal,true);
  run_simulation();
  // collect samples
  if (mpi_comm.is_master()) {
    for (const mpi::proc& p : working_procs) {
      //std::cout << "reciving results from p = " << p << "\n";
      MPI_recv_results(mpi_comm, p, mpi::MP_data_samples);
    }
    finalize_results();
  	print_results();
  }
  else {
    //std::cout << "sending results from p = " << mpi_comm.rank() << "\n";
    MPI_send_results(mpi_comm, mpi_comm.master(), mpi::MP_data_samples);
  }
  return 0;
}

int VMCRun::master_run(const Eigen::VectorXd& vparms, double& en_mean, double& en_stddev,
  Eigen::VectorXd& grad, Eigen::MatrixXd& sr_matrix, const mpi::mpi_communicator& mpi_comm)
{
  bool with_psi_grad = true;
  build_config(vparms, with_psi_grad);
  reset_observables();
  do_warmup_run();
  do_measure_run(num_samples_);
  // recieve from works
  for (const mpi::proc& p : working_procs) {
    MPI_recv_results(mpi_comm, p, mpi::MP_data_samples);
  }
  finalize_results();
}

int VMCRun::worker_run(const mpi::mpi_communicator& mpi_comm)
{
  while (true) {
  	mpi::mpi_status msg = mpi_comm.probe();
  	switch (msg:tag()) {
  	  case mpi::MP_start_simulation:
  		bool with_psi_grad = true;
  		build_config(vparms, with_psi_grad);
  		reset_observables();
  		do_warmup_run();
  		do_measure_run(num_samples_);
    	MPI_send_results(mpi_comm, mpi_comm.master(), mpi::MP_data_samples);
  	  	break;
  	  case mpi::MP_stop_simulation:
        mpi_comm.recv(msg.source(), msg.tag());
        return 0;
  	  default: break;
  	}
  }
  return 0;
}
#endif


} // end namespace vmc
/*---------------------------------------------------------------------------
* @Author: Amal Medhi
* @Date:   2022-10-15 14:57:27
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2022-10-15 14:57:27
* Copyright (C) 2015-2022 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
*----------------------------------------------------------------------------*/
#ifndef VMCRUN_H
#define VMCRUN_H

#include <iostream>
#include <set>
#include "./vmc.h"

namespace vmc {

enum class mpi_mode {NORMAL, BC_TWIST};

class VMCRun : public VMC
{
public:
  VMCRun(const input::Parameters& parms); 
  ~VMCRun() {}
  int run(const input::Parameters& inputs); 
  int run(const input::Parameters& inputs, const mpi::mpi_communicator& mpi_comm);
  int run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err, RealMatrix& sr_matrix, const bool& silent_mode=false);
  int run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err, const bool& silent_mode=false);
  // external optimizer
  //bool Evaluate(const double* parameters, double* cost, double* gradient) override;
  //int NumParameters() const override;

  int start_worker_run(void);
  int stop_worker_run(void);
  //int start_worker_run(const mpi::mpi_communicator& mpi_comm);
  void set_mpi_comm(const mpi::mpi_communicator& mpi_comm) { mpi_comm_=mpi_comm; }
  void set_measure_steps(const int& n) { measure_samples_ = n; }
  void set_workers_list(const std::vector<mpi::proc>& workers) { worker_procs_=workers; }
  void set_bc_list(const std::vector<int>& bc_list) { bc_list_=bc_list; }
  void set_mpi_mode(const mpi_mode& mode) { mpi_mode_=mode; }
private:
  // for MPI run
  mpi::mpi_communicator mpi_comm_;
  mpi_mode mpi_mode_{mpi_mode::NORMAL};
  std::vector<mpi::proc> worker_procs_; // only for master
  std::vector<int> bc_list_;
  int measure_samples_{0};

  int master_run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err, RealMatrix& sr_matrix, const bool& silent_mode=false);
  int master_run(const var::parm_vector& vparms, double& en, double& en_err,
	RealVector& grad, RealVector& grad_err, const bool& silent_mode=false);
};


} // end namespace vmc

#endif
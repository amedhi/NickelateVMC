/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <iostream>
#include "./vmc.h"
//#include "./stochastic_reconf.h"
#include "./optimizer.h"
//#include "../optimizer/optimizer.h"

namespace vmc {

class Simulator : public scheduler::Worker
{
public:
  Simulator(const input::Parameters& parms); 
  ~Simulator() {}
  int start(const mpi::mpi_communicator& mpi_comm) override; 
  int run(const input::Parameters& parms) override;
  int run(const input::Parameters& parms, 
    const mpi::mpi_communicator& mpi_comm) override;
  void finish(void) override {} 
  void dostep(void) override {} 
  void halt(void) override {} 
  static void print_copyright(std::ostream& os);
  //const VMC& sim(void) { return vmc; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  VMCRun vmc_;
  Optimizer optimizer_;
  //StochasticReconf sreconf;
  //opt::Optimizer optimizer_;
  mpi_mode mpi_mode_; 
  bool optimization_mode_{false};
  int mpirun_vmc(const mpi::mpi_communicator& mpi_comm, 
    const std::set<mpi::proc>& working_procs, const int& proc_samples,
    const bool& quiet=false); 
  int mpirun_vmc(const mpi::mpi_communicator& mpi_comm, 
    const std::set<mpi::proc>& working_procs, const std::vector<int>& bc_list,
    const bool& quiet=false); 
  //static double enfunc(const std::vector<double>& x, std::vector<double>& grad, 
  //  void *my_func_data);
};

//double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data);

} // end namespace vmc

#endif

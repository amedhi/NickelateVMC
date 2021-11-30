/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:14:38
*----------------------------------------------------------------------------*/
#include <iomanip>
#include <set>
#include <boost/algorithm/string.hpp>
#include "simulator.h"
//#include <nlopt.hpp>
//#include "../optimizer/LBFGS.h"

namespace vmc {

Simulator::Simulator(const input::Parameters& inputs) : vmc(inputs)
{
  optimization_mode_ = inputs.set_value("optimizing_run",false);
  if (optimization_mode_) {
    sreconf.init(inputs, vmc);
    //optimizer_.init(inputs, vmc);
    //vmc.set_box_constraints();
    //nlopt_.init(inputs, vmc);
  }
  // mpi mode
  int nowarn;
  std::string mode = inputs.set_value("mpi_mode","NORMAL",nowarn);
  boost::to_upper(mode);
  if (mode=="BC_TWIST") {
    mpi_mode_ = mpi_mode::BC_TWIST; 
  }
  else {
    mpi_mode_ = mpi_mode::NORMAL; 
  }
}

int Simulator::run(const input::Parameters& inputs) 
{
  // optimization run
  if (optimization_mode_) {
    if (!inputs.have_option_quiet()) std::cout << " starting optimizing run\n";
    //vmc.start(inputs, run_mode::energy_function, false);
    //nlopt_.optimize(vmc);
    vmc.start(inputs, run_mode::sr_function, true);
    if (sreconf.optimize(vmc)) {
    //if (optimizer_.optimize(vmc)) {
      vmc.run_simulation(sreconf.optimal_parms());
      vmc.print_results();
    }
    return 0;
  }

  // normal run
  if (!inputs.have_option_quiet()) std::cout << " starting vmc run\n";
  vmc.start(inputs, run_mode::normal);
  vmc.run_simulation();
  vmc.print_results();
  return 0;
}

// parallel run
int Simulator::run(const input::Parameters& inputs, 
  const mpi::mpi_communicator& mpi_comm)
{
  int num_procs = mpi_comm.size();
  std::set<mpi::proc> working_procs;  
  /**************************************************
   * NORMAL MODE: Split samples into different procs
   **************************************************/
  if (mpi_mode_ == mpi_mode::NORMAL) {
    int total_samples = vmc.num_measure_steps();
    int proc_samples = 0;
    if (total_samples >= num_procs) {
      // distribute equally (master takes lesser load)
      proc_samples = total_samples/num_procs;
      int excess_samples = total_samples-proc_samples*num_procs;
      if (mpi_comm.is_master()) {
        proc_samples -= num_procs-excess_samples-1;
      }
      else {
        proc_samples += 1;
      }
      // save the working slaves
      if (mpi_comm.is_master()) {
        for (const mpi::proc& p : mpi_comm.slave_procs()) {
          working_procs.insert(p);
        }
      }
    }
    else {
      if (mpi_comm.is_master()) {
        std::cout << ">> alert: redundant processes in MPI run\n"; 
      }
      if (mpi_comm.rank() >= total_samples) return 0; // no work for you
      proc_samples = 1; // for the remaining guys
      // save the working slaves
      if (mpi_comm.is_master()) {
        int i = 0;
        for (const mpi::proc& p : mpi_comm.slave_procs()) {
          if (++i == total_samples) break;
          working_procs.insert(p);
        }
      }
    }
    // do one MPI run
    if (mpi_comm.is_master()) {
      if (!inputs.have_option_quiet()) std::cout << " starting vmc run in normal MPI mode\n";
    }
    mpirun_vmc(inputs, mpi_comm, working_procs, proc_samples, 
      run_mode::normal);
    if (mpi_comm.is_master()) {
      vmc.print_results();
    }
    return 0;
  }

  /*******************************************************
   * BC_TWISTS: Run different BC_TWIST in different procs
   *******************************************************/
  if (mpi_mode_ == mpi_mode::BC_TWIST) {
    int num_bcs = vmc.num_boundary_twists();
    std::vector<int> bc_list;
    if (num_bcs >= num_procs) {
      int num_proc_bcs = num_bcs/num_procs;
      if (mpi_comm.is_master()) {
        for (int i=0; i<num_proc_bcs; ++i) bc_list.push_back(i);
      }
      else {
        if (num_proc_bcs*num_procs != num_bcs) num_proc_bcs++;
        for (int i=0; i<num_proc_bcs; ++i) {
          int id = mpi_comm.rank()*num_proc_bcs+i;
          if (id >= num_bcs) break;
          bc_list.push_back(id);
        }
      }
      // working procs
      if (mpi_comm.is_master()) {
        for (const mpi::proc& p : mpi_comm.slave_procs()) {
          working_procs.insert(p);
        }
      }
    }
    else {
      if (mpi_comm.is_master()) {
        std::cout << ">> alert: redundant processes in MPI run\n"; 
      }
      if (mpi_comm.rank() >= num_bcs) return 0; // no work for you
      bc_list.push_back(mpi_comm.rank());
      // save the working slaves
      if (mpi_comm.is_master()) {
        int i = 0;
        for (const mpi::proc& p : mpi_comm.slave_procs()) {
          if (++i == num_bcs) break;
          working_procs.insert(p);
        }
      }
    }
    // do one MPI run
    if (mpi_comm.is_master()) {
      if (!inputs.have_option_quiet()) std::cout << " starting vmc run in BC_TWIST MPI mode\n";
    }
    mpirun_vmc(inputs, mpi_comm, working_procs, bc_list, run_mode::normal);
    if (mpi_comm.is_master()) {
      vmc.print_results();
    }
    return 0;
    //vmc.start(inputs, run_mode::normal);
    //vmc.run_simulation();
    //if (mpi_comm.is_master()) {
    //  vmc.print_results();
    //}
  }

  return 0;


  // disordered system
  if (vmc.disordered_system()) {
    int num_proc = mpi_comm.size();
    int num_dconf = vmc.num_disorder_configs();
    int n1, n2;
    if (num_proc==num_dconf) {
      n1 = mpi_comm.rank();
      n2 = mpi_comm.rank();
    }
    else if (num_proc > num_dconf) {
      n1 = mpi_comm.rank();
      n2 = mpi_comm.rank();
      if (n1 >= num_dconf) return 0; // no job for you
    }
    else {
      int n = std::nearbyint(static_cast<double>(num_dconf)/num_proc);
      n1 = n * mpi_comm.rank();
      n2 = n1 + n;
      if (n1 >= num_dconf) return 0; // no job for you
      if (mpi_comm.rank()==(num_proc-1)) {
        n2 = num_dconf;
      }
    }
    // optimizing run
    if (optimization_mode_) {
      for (unsigned n=n1; n<n2; ++n) {
        if (vmc.optimal_parms_exists(n)) continue;
        std::cout << " optimizing disorder config " << n;
        std::cout << " of " << vmc.num_disorder_configs() << "\n";
        vmc.start(inputs, run_mode::sr_function, true);
        vmc.set_disorder_config(n);
        if (sreconf.optimize(vmc)) {
          vmc.save_optimal_parms(sreconf.optimal_parms());
          //vmc.run_simulation(sreconf.optimal_parms());
          //vmc.print_results();
        }
      }
    }
  }
  return 0;
}

// one run in parallel
int Simulator::mpirun_vmc(const input::Parameters& inputs, 
  const mpi::mpi_communicator& mpi_comm, 
  const std::set<mpi::proc>& working_procs,
  const int& proc_samples, const run_mode& mode) 
{
  //std::cout << "starting in p = " << mpi_comm.rank() << "\n";
  vmc.start(inputs, mode);
  //std::cout << "warming in p = " << mpi_comm.rank() << "\n";
  vmc.do_warmup_run();
  //std::cout << "measuring "<< proc_samples<<" in p = " << mpi_comm.rank() << "\n";
  vmc.do_measure_run(proc_samples);
  // collect samples
  if (mpi_comm.is_master()) {
    /*
    while (true) {
      mpi::mpi_status msg = mpi_comm.probe();
      switch (msg.tag()) {
        case mpi::MP_data_samples:
          vmc.MPI_recv_results(mpi_comm, msg.source(), msg.tag());
          break;
        default:
          std::cout << "Simulator::run: unexpected message\n";
          return 1;
      }
    }*/
    for (const mpi::proc& p : working_procs) {
      //std::cout << "reciving results from p = " << p << "\n";
      vmc.MPI_recv_results(mpi_comm, p, mpi::MP_data_samples);
    }
    vmc.finalize_results();
  }
  else {
    //std::cout << "sending results from p = " << mpi_comm.rank() << "\n";
    vmc.MPI_send_results(mpi_comm, mpi_comm.master(), mpi::MP_data_samples);
  }
  return 0;
}

int Simulator::mpirun_vmc(const input::Parameters& inputs, const mpi::mpi_communicator& 
  mpi_comm, const std::set<mpi::proc>& working_procs, 
  const std::vector<int>& bc_list, const run_mode& mode) 
{
  vmc.start(inputs, mode, true);
  if (!inputs.have_option_quiet()) {
    std::cout << " running BC_TWISTS: "<<bc_list.front()+1<<"-"
      <<bc_list.back()+1<<" in p"<<mpi_comm.rank()<<"\n";
  } 
  vmc.run_simulation(-1, bc_list);
  if (mpi_comm.is_master()) {
    //if (!inputs.have_option_quiet()) {
    //  std::cout << " accumulating results...\n";
    //}
    for (const mpi::proc& p : working_procs) {
      vmc.MPI_recv_results(mpi_comm, p, mpi::MP_data_samples);
    }
  }
  else {
    vmc.MPI_send_results(mpi_comm, mpi_comm.master(), mpi::MP_data_samples);
  }
  return 0;
}



  /*
  //nlopt::opt opt(nlopt::LD_MMA, varparms.size());
  nlopt::opt opt(nlopt::LD_LBFGS, varparms.size());
  opt.set_lower_bounds(simulator.vparm_lbound());
  opt.set_upper_bounds(simulator.vparm_ubound());
  opt.set_xtol_rel(1e-3);
  opt.set_min_objective(wrapper, &simulator);
  double mine;
  nlopt::result result = opt.optimize(varparms, mine);
  //std::cout << nlopt::result << "\n";
  */

  /*using namespace LBFGSpp;
  LBFGSParam<double> lbfgs_param;
  LBFGSSolver<double> solver(lbfgs_param);
  //simulator.run(varparms, true);
  double emin;
  int niter = solver.minimize(simulator, varparms, emin);
  std::cout << niter << " iterations" << std::endl;
  std::cout << "x = \n" << varparms.transpose() << std::endl;
  std::cout << "f(x) = " << emin << std::endl;
  */
  
  /*
  std::cout << "var parms = " << varparms_.size();
  simulator.optimizing(varparms_,energy,energy_grad);
  int j = 0;
  for (int i=0; i<10; ++i) {
    simulator.run(varparms_,energy,energy_grad);
    varparms_[j++] += 0.2;
    if (j == varparms_.size()) j = 0;
  }
  return 0;
  */

/*
double wrapper(const std::vector<double>& x, std::vector<double>& grad, void *my_data)
{
  Eigen::VectorXd en_grad(x.size());
  double en = reinterpret_cast<Simulator*>(my_data)->energy_function(x,en_grad);
  for (int i=0; i<grad.size(); ++i) grad[i] = en_grad[i];
  return en;
}
*/

void Simulator::print_copyright(std::ostream& os)
{
  VMC::copyright_msg(os);
}


} // end namespace mc

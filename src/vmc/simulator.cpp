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

Simulator::Simulator(const input::Parameters& inputs) : vmc_(inputs)
{
  optimization_mode_ = inputs.set_value("optimizing_run",false);
  if (optimization_mode_) {
    optimizer_.init(inputs, vmc_);
  }
}

int Simulator::start(const mpi::mpi_communicator& mpi_comm) 
{
  if (optimization_mode_) {
    vmc_.set_mpi_comm(mpi_comm);
    if (mpi_comm.is_master()) optimizer_.print_info(vmc_);
  }

  // Job division in case of MPI run
  int num_procs = mpi_comm.size();
  if (num_procs > 1) {
    // MPI Mode
    if (vmc_.num_boundary_twists()>1) {
      mpi_mode_ = mpi_mode::BC_TWIST; 
    }
    else {
      mpi_mode_ = mpi_mode::NORMAL; 
    }
    vmc_.set_mpi_mode(mpi_mode_);

    /**************************************************
    * NORMAL MODE: Split samples into different procs
    **************************************************/
    if (mpi_mode_ == mpi_mode::NORMAL) {
      std::vector<mpi::proc> worker_list;  
      int total_samples = vmc_.num_measure_steps();
      int proc_samples = 0;
      if (total_samples >= num_procs) {
        // distribute equally (master takes lesser load)
        proc_samples = total_samples/num_procs;
        int excess_samples = total_samples-proc_samples*num_procs;
        // assign the excess samples (except for the master)
        if (!mpi_comm.is_master() && mpi_comm.rank()<=excess_samples) {
          proc_samples += 1;
        }
        /*
        if (mpi_comm.is_master()) {
          proc_samples -= num_procs-excess_samples-1;
        }
        else {
          proc_samples += 1;
        }*/
        vmc_.set_measure_steps(proc_samples);
        // list of workers 
        if (mpi_comm.is_master()) {
          for (const mpi::proc& p : mpi_comm.slave_procs()) {
            worker_list.push_back(p);
          }
          vmc_.set_workers_list(worker_list);
        }
      }
      else {
        if (mpi_comm.is_master()) {
          std::cout << " >> alert: redundant processes in MPI run\n"; 
        }
        if (mpi_comm.rank() >= total_samples) return 0; // no work for you
        proc_samples = 1; // for the remaining guys
        vmc_.set_measure_steps(proc_samples);
        // list of workers 
        if (mpi_comm.is_master()) {
          int i = 0;
          for (const mpi::proc& p : mpi_comm.slave_procs()) {
            if (++i == total_samples) break;
            worker_list.push_back(p);
          }
          vmc_.set_workers_list(worker_list);
        }
      }
    }

    /*******************************************************
     * BC_TWISTS: Run different BC_TWIST in different procs
     *******************************************************/
    if (mpi_mode_ == mpi_mode::BC_TWIST) {
      int num_bcs = vmc_.num_boundary_twists();
      std::vector<int> bc_list;
      std::vector<mpi::proc> worker_list;  

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
            worker_list.push_back(p);
          }
          vmc_.set_workers_list(worker_list);
        }
      }
      else {
        if (mpi_comm.is_master()) {
          std::cout << " >> alert: redundant processes in MPI run\n"; 
        }
        if (mpi_comm.rank() >= num_bcs) return 0; // no work for you
        bc_list.push_back(mpi_comm.rank());
        // save the working slaves
        if (mpi_comm.is_master()) {
          int i = 0;
          for (const mpi::proc& p : mpi_comm.slave_procs()) {
            if (++i == num_bcs) break;
            worker_list.push_back(p);
          }
          vmc_.set_workers_list(worker_list);
        }
      }
      // set the list of BCs to run
      vmc_.set_bc_list(bc_list);
    } // end BC_TWIST
  } // end if mpirun
  return 0;
}

// serial run
int Simulator::run(const input::Parameters& inputs) 
{
  if (!optimization_mode_) {
    return vmc_.run(inputs);
  }

  // Optimization
  optimizer_.optimize(vmc_);
  return 0;
}

// parallel run
int Simulator::run(const input::Parameters& inputs, 
  const mpi::mpi_communicator& mpi_comm)
{
  if (!optimization_mode_) {
    return vmc_.run(inputs,mpi_comm);
  }

  // Optimization
  if (mpi_comm.is_master()) {
    if (!inputs.have_option_quiet()) std::cout << " starting optimizing run\n";
    optimizer_.optimize(vmc_);
    // terminate workers
    vmc_.stop_worker_run();
  }
  else {
    vmc_.start_worker_run();
  }
  // all should meet here, before the next run
  mpi_comm.barrier(); 
  return 0;
}

//-------------------------
#ifdef ON
//-------------------------


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

    // Optimizing run
    if (optimization_mode_) {
      if (mpi_comm.is_master()) {
        if (!inputs.have_option_quiet()) std::cout << " starting optimizing run\n";
      }
      int num_vparms = vmc.num_varp();
      Eigen::VectorXd grad(num_vparms);
      Eigen::MatrixXd sr_matrix(num_vparms,num_vparms);
      // start vmc in sr_function mode in all procs
      vmc.start(inputs, run_mode::sr_function, true);
      for (int n=0; n<sreconf.num_opt_samples(); ++n) {
        // starting parameters for each sample
        var::parm_vector vparms = vmc.varp_values();
        if (mpi_comm.is_master()) sreconf.start(vmc, n);
        while (true) {
          mpi_comm.barrier(); // for safety
          bool with_psi_grad = true;
          vmc.build_config(vparms, with_psi_grad);
          mpirun_vmc(mpi_comm,working_procs,proc_samples,true); 
          if (mpi_comm.is_master()) {
            // quantities needed for optimization
            double energy = vmc.observable().energy().mean();
            double energy_err = vmc.observable().energy().stddev();
            grad = vmc.observable().energy_grad().mean_data();
            vmc.observable().sr_matrix().get_matrix(sr_matrix);
            // iterate
            exit_stat status;
            bool done=sreconf.iterate(vparms,energy,energy_err,grad,sr_matrix,
              status);
            if (done) {
              sreconf.finalize();
              for (const mpi::proc& p : working_procs) {
                mpi_comm.send(p, mpi::MP_stop_simulation);
              }
              break;
            }
            else {
              for (const mpi::proc& p : working_procs) {
                mpi_comm.send(p, mpi::MP_variational_parms, vparms);
              }
            }
          }
          else {
            mpi::mpi_status msg = mpi_comm.probe();
            if (msg.tag() == mpi::MP_variational_parms) {
              mpi_comm.recv(msg.source(), msg.tag(), vparms);
            }
            else if (msg.tag() == mpi::MP_stop_simulation) {
              mpi_comm.recv(msg.source(), msg.tag());
              break;
            }
            else {
              std::cout << "Simulator::run: unexpected message in optimizing run\n";
            }
          }
        } // while each sample
      }// all samples
    }// optimizing run

    // VMC run
    else {
      if (mpi_comm.is_master()) {
        if (!inputs.have_option_quiet()) std::cout << " starting vmc run in normal MPI mode\n";
      }
      vmc.start(inputs,run_mode::normal,true);
      mpirun_vmc(mpi_comm,working_procs,proc_samples); 
      if (mpi_comm.is_master()) {
        vmc.print_results();
      }
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

    // Optimizing run
    if (optimization_mode_) {
      if (mpi_comm.is_master()) {
        if (!inputs.have_option_quiet()) std::cout << " starting optimizing run\n";
      }
      int num_vparms = vmc.num_varp();
      Eigen::VectorXd grad(num_vparms);
      Eigen::MatrixXd sr_matrix(num_vparms,num_vparms);
      // start vmc in sr_function mode in all procs
      vmc.start(inputs,run_mode::sr_function,true);
      for (int n=0; n<sreconf.num_opt_samples(); ++n) {
        // starting parameters for each sample
        var::parm_vector vparms = vmc.varp_values();
        if (mpi_comm.is_master()) sreconf.start(vmc, n);
        while (true) {
          mpi_comm.barrier(); // for safety
          bool with_psi_grad = true;
          vmc.build_config(vparms, with_psi_grad);
          mpirun_vmc(mpi_comm,working_procs,bc_list,true);
          if (mpi_comm.is_master()) {
            // quantities needed for optimization
            double energy = vmc.observable().energy().mean();
            double energy_err = vmc.observable().energy().stddev();
            grad = vmc.observable().energy_grad().mean_data();
            vmc.observable().sr_matrix().get_matrix(sr_matrix);
            // iterate
            exit_stat status;
            bool done=sreconf.iterate(vparms,energy,energy_err,grad,sr_matrix,
              status);
            if (done) {
              sreconf.finalize();
              for (const mpi::proc& p : working_procs) {
                mpi_comm.send(p, mpi::MP_stop_simulation);
              }
              break;
            }
            else {
              for (const mpi::proc& p : working_procs) {
                mpi_comm.send(p, mpi::MP_variational_parms, vparms);
              }
            }
          }
          else {
            mpi::mpi_status msg = mpi_comm.probe();
            if (msg.tag() == mpi::MP_variational_parms) {
              mpi_comm.recv(msg.source(), msg.tag(), vparms);
            }
            else if (msg.tag() == mpi::MP_stop_simulation) {
              mpi_comm.recv(msg.source(), msg.tag());
              break;
            }
            else {
              std::cout << "Simulator::run: unexpected message in optimizing run\n";
            }
          }
        } // while each sample
      }// all samples
    }// optimizing run

    // VMC run
    else {
      if (mpi_comm.is_master()) {
        if (!inputs.have_option_quiet()) std::cout << " starting vmc run in BC_TWIST MPI mode\n";
      }
      vmc.start(inputs,run_mode::normal,true);
      mpirun_vmc(mpi_comm,working_procs,bc_list,inputs.have_option_quiet());
      if (mpi_comm.is_master()) {
        vmc.print_results();
      }
      return 0;
    }
  }
  return 0;
}

// one run in parallel
int Simulator::mpirun_vmc(const mpi::mpi_communicator& mpi_comm, 
  const std::set<mpi::proc>& working_procs, const int& proc_samples,
  const bool& quiet)
{
  if (!quiet) {
    std::cout << " running "<<proc_samples<<" samples in p"<<mpi_comm.rank()<<"\n";
  } 

  vmc.reset_observables();
  //std::cout << "warming up in p" << mpi_comm.rank() << "\n";
  vmc.do_warmup_run();
  //std::cout << "measuring "<<proc_samples<<" in p" << mpi_comm.rank() << "\n";
  vmc.do_measure_run(proc_samples);
  // collect samples
  if (mpi_comm.is_master()) {
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

int Simulator::mpirun_vmc(const mpi::mpi_communicator& mpi_comm, 
  const std::set<mpi::proc>& working_procs, const std::vector<int>& bc_list,
  const bool& quiet) 
{
  if (!quiet) {
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

//-------------------------
#endif
//-------------------------


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

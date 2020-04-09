/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-09 15:27:50
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:14:38
*----------------------------------------------------------------------------*/
#include <iomanip>
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
}

int Simulator::run(const input::Parameters& inputs) 
{
  // disordered system
  if (vmc.disordered_system()) {
    // optimizing run
    if (optimization_mode_) {
      for (int n=0; n<vmc.num_disorder_configs(); ++n) {
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
      return 0;
    }

    // normal run
    for (int n=0; n<vmc.num_disorder_configs(); ++n) {
    //for (unsigned n=0; n<1; ++n) {
      vmc.disorder_start(inputs, n);
      vmc.run_simulation();
    //  vmc.start(inputs, n);
    //  vmc.run_simulation();
      //vmc.save_results();
    }
    //vmc.print_avg_results();
    vmc.print_results();
    return 0;
  }

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
  const scheduler::mpi_communicator& mpi_comm)
{
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

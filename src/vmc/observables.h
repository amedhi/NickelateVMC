/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#ifndef MC_OBSERVABLES_H
#define MC_OBSERVABLES_H

#include <string>
#include <vector>
#include <stdexcept>
#include "../scheduler/task.h"
#include "../lattice/lattice.h"
#include "../lattice/graph.h"
#include "../model/model.h"
#include "./sysconfig.h"
#include "./energy.h"
#include "./sccorr.h"
#include "./particle.h"

namespace vmc {

class ObservableSet 
{
public:
  ObservableSet();
  ~ObservableSet() {}
  std::stringstream& headstream(void) { return headstream_; }
  //void init(const input::Parameters& inputs, 
  //  void (&print_copyright)(std::ostream& os), const lattice::LatticeGraph& graph, 
  //  const model::Hamiltonian& model, const SysConfig& config);
  void init(const input::Parameters& inputs, const lattice::LatticeGraph& graph, 
    const model::Hamiltonian& model, const SysConfig& config, const std::string& prefix);
  void as_functions_of(const std::vector<std::string>& xvars=std::vector<std::string>());
  void as_functions_of(const std::string& xvar);
  void switch_off(void);
  void reset(void); 
  void reset_grand_data(void); 
  void save_results(void); 
  void avg_grand_data(void); 
  int do_measurement(const lattice::LatticeGraph& graph, 
    const model::Hamiltonian& model, const SysConfig& config, const SiteDisorder& site_disorder);
  inline Energy& energy(void) { return energy_; }
  inline EnergyGradient& energy_grad(void) { return energy_grad_; }
  inline SC_Correlation& sc_corr(void) { return sc_corr_; }
  inline SR_Matrix& sr_matrix(void) { return sr_matrix_; }
  void finalize(void);
  void print_heading(void);
  void print_results(const std::vector<double>& xvals=std::vector<double>()); 
  void print_results(const double& xval); 

  void MPI_send_results(const mpi::mpi_communicator& mpi_comm, 
    const mpi::proc& proc, const int& msg_tag);
  void MPI_recv_results(const mpi::mpi_communicator& mpi_comm, 
    const mpi::proc& proc, const int& msg_tag); 
private:
  bool replace_mode_{true};
  std::stringstream headstream_;
  std::vector<std::string> xvars_;
  unsigned num_xvars_{0};
  Energy energy_;
  EnergyGradient energy_grad_;
  SC_Correlation sc_corr_;
  SR_Matrix sr_matrix_;
  SiteOccupancy site_occupancy_;
};





} // end namespace vmc

#endif

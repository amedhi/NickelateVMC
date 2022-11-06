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
#include "../model/model.h"
#include "./sysconfig.h"
#include "./energy.h"
#include "./particle.h"
#include "./spincorr.h"
#include "./sccorr.h"

namespace vmc {

class ObservableSet 
{
public:
  ObservableSet();
  ~ObservableSet() {}
  std::stringstream& headstream(void) { return headstream_; }
  //void init(const input::Parameters& inputs, 
  //  void (&print_copyright)(std::ostream& os), const lattice::Lattice& lattice, 
  //  const model::Hamiltonian& model, const SysConfig& config);
  void init(const input::Parameters& inputs, const lattice::Lattice& lattice, 
    const model::Hamiltonian& model, const SysConfig& config, const std::string& prefix,
    const int& sample_size);
  void as_functions_of(const std::vector<std::string>& xvars=std::vector<std::string>());
  void as_functions_of(const std::string& xvar);
  void switch_off(void);
  void reset(void); 
  void reset_batch_limit(const int& sample_size);
  void reset_grand_data(void); 
  void save_results(void); 
  void avg_grand_data(void); 
  int do_measurement(const lattice::Lattice& lattice, 
    const model::Hamiltonian& model, const SysConfig& config, const SiteDisorder& site_disorder);
  inline Energy& energy(void) { return energy_; }
  inline EnergyGradient& energy_grad(void) { return energy_grad_; }
  inline SC_Correlation& sc_corr(void) { return sc_corr_; }
  inline SR_Matrix& sr_matrix(void) { return sr_matrix_; }
  const Energy& energy(void) const { return energy_; }
  const EnergyGradient& energy_grad(void) const { return energy_grad_; }
  const SpinCorrelation& spin_corr(void) const { return spin_corr_; }
  const SC_Correlation& sc_corr(void) const { return sc_corr_; }
  const SR_Matrix& sr_matrix(void) const { return sr_matrix_; }
  const ParticleDensity& particle_density(void) const { return particle_density_; }

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
  ParticleDensity particle_density_;
  SpinCorrelation spin_corr_;
  SC_Correlation sc_corr_;
  SR_Matrix sr_matrix_;
};





} // end namespace vmc

#endif

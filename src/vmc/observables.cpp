/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
*----------------------------------------------------------------------------*/
#include "observables.h"
#include <boost/algorithm/string.hpp>

namespace vmc {

ObservableSet::ObservableSet() 
  : energy_("Energy")
  , energy_grad_("EnergyGradient")
  , spin_corr_("SpinCorrelation")
  , sc_corr_("SC_Correlation")
  , sr_matrix_("SR_Matrix")
  , particle_density_("ParticleDensity")
  , doublon_density_("DoublonDensity")
{
}

//void ObservableSet::init(const input::Parameters& inputs, 
//    void (&print_copyright)(std::ostream& os), const lattice::Lattice& lattice, 
//    const model::Hamiltonian& model, const SysConfig& config)
void ObservableSet::init(const input::Parameters& inputs, 
  const lattice::Lattice& lattice, const model::Hamiltonian& model, 
  const SysConfig& config, const std::string& prefix, const int& sample_size)
{
  // file open mode
  std::string mode = inputs.set_value("mode", "NEW");
  boost::to_upper(mode);
  if (mode=="APPEND") replace_mode_ = false;
  else replace_mode_ = true;
  // check which observables to calculate
  //for (auto& obs : *this) obs.get().check_on(inputs, replace_mode_);
  // heading message
  //print_copyright(headstream_);
  //model.print_info(headstream_);
  num_xvars_ = 0; 

  // files
  energy_.set_ofstream(prefix);
  energy_grad_.set_ofstream(prefix);
  spin_corr_.set_ofstream(prefix);
  sc_corr_.set_ofstream(prefix);
  sr_matrix_.set_ofstream(prefix);
  particle_density_.set_ofstream(prefix);
  doublon_density_.set_ofstream(prefix);

  // switch on required observables
  energy_.check_on(inputs, replace_mode_);
  energy_grad_.check_on(inputs, replace_mode_);
  if (energy_grad_) energy_.switch_on();
  spin_corr_.check_on(inputs,replace_mode_);
  sc_corr_.check_on(inputs,replace_mode_);
  sr_matrix_.check_on(inputs,replace_mode_);
  particle_density_.check_on(inputs,replace_mode_);
  doublon_density_.check_on(inputs,replace_mode_);

  // set up observables
  if (energy_) energy_.setup(lattice,model);
  if (energy_grad_) energy_grad_.setup(config,sample_size);
  if (spin_corr_) spin_corr_.setup(lattice);
  if (sc_corr_) sc_corr_.setup(lattice,config.wavefunc().mf_order(),sample_size);
  if (sr_matrix_) sr_matrix_.setup(lattice,config);
  if (particle_density_) particle_density_.setup(lattice,config);
  if (doublon_density_) doublon_density_.setup(lattice,config);
}

void ObservableSet::reset(void)
{
  if (energy_) energy_.reset();
  if (energy_grad_) energy_grad_.reset();
  if (spin_corr_) spin_corr_.reset();
  if (sc_corr_) sc_corr_.reset();
  if (sr_matrix_) sr_matrix_.reset();
  if (particle_density_) particle_density_.reset();
  if (doublon_density_) doublon_density_.reset();
}

void ObservableSet::reset_batch_limit(const int& sample_size)
{
  if (energy_grad_) energy_grad_.reset_batch_limit(sample_size);
  if (sc_corr_) sc_corr_.reset_batch_limit(sample_size);
  if (sr_matrix_) sr_matrix_.reset_batch_limit(sample_size);
}

void ObservableSet::reset_grand_data(void)
{
  if (energy_) energy_.reset_grand_data();
  if (energy_grad_) energy_grad_.reset_grand_data();
  if (spin_corr_) spin_corr_.reset_grand_data();
  if (sc_corr_) sc_corr_.reset_grand_data();
  if (sr_matrix_) sr_matrix_.reset_grand_data();
  if (particle_density_) particle_density_.reset_grand_data();
  if (doublon_density_) doublon_density_.reset_grand_data();
}

void ObservableSet::save_results(void)
{
  if (energy_) energy_.save_result();
  if (energy_grad_) energy_grad_.save_result();
  if (spin_corr_) spin_corr_.save_result();
  if (sc_corr_) sc_corr_.save_result();
  if (sr_matrix_) sr_matrix_.save_result();
  if (doublon_density_) doublon_density_.save_result();
}

void ObservableSet::avg_grand_data(void)
{
  if (energy_) energy_.avg_grand_data();
  if (energy_grad_) energy_grad_.avg_grand_data();
  if (spin_corr_) spin_corr_.avg_grand_data();
  if (sc_corr_) sc_corr_.avg_grand_data();
  if (sr_matrix_) sr_matrix_.avg_grand_data();
  if (particle_density_) particle_density_.avg_grand_data();
  if (doublon_density_) doublon_density_.avg_grand_data();
}


int ObservableSet::do_measurement(const lattice::Lattice& lattice, 
    const model::Hamiltonian& model, const SysConfig& config, const SiteDisorder& site_disorder)
{
  if (energy_) energy_.measure(lattice,model,config,site_disorder);
  if (energy_grad_) {
    if (!energy_) 
      throw std::logic_error("ObservableSet::measure: dependency not met for 'energy'");
    energy_grad_.measure(config, energy_.config_value().sum());
  }
  if (spin_corr_) spin_corr_.measure(lattice,model,config);
  if (sc_corr_) sc_corr_.measure(lattice,model,config);
  if (sr_matrix_) {
    if (!energy_grad_) 
      throw std::logic_error("ObservableSet::measure: dependency not met for 'sr_matrix_'");
    sr_matrix_.measure(energy_grad_.grad_logpsi());
  }
  if (particle_density_) particle_density_.measure(lattice, config);
  if (doublon_density_) doublon_density_.measure(lattice, config);
  return 0;
}

void ObservableSet::finalize(void)
{
  if (energy_grad_) {
    energy_grad_.finalize();
  }
}

void ObservableSet::as_functions_of(const std::vector<std::string>& xvars)
{
  xvars_ = xvars;
  num_xvars_ = xvars_.size();
}

void ObservableSet::as_functions_of(const std::string& xvar)
{
  xvars_ = {xvar};
  num_xvars_ = 1;
}

void ObservableSet::switch_off(void) {
  energy_.switch_off();
  energy_grad_.switch_off();
  spin_corr_.switch_off();
  sc_corr_.switch_off();
  sr_matrix_.switch_off();
  particle_density_.switch_off();
  doublon_density_.switch_off();
}

void ObservableSet::print_heading(void)
{
  energy_.print_heading(headstream_.rdbuf()->str(),xvars_);
  energy_grad_.print_heading(headstream_.rdbuf()->str(),xvars_);
  spin_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
  sc_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
  sr_matrix_.print_heading(headstream_.rdbuf()->str(),xvars_);
  particle_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
  doublon_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
}

void ObservableSet::print_results(const std::vector<double>& xvals) 
{
  if (num_xvars_ != xvals.size()) 
    throw std::invalid_argument("Observables::print_result: 'x-vars' size mismatch");
  if (energy_) {
    energy_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_.print_result(xvals);
  }
  if (energy_grad_) {
    energy_grad_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_grad_.print_result(xvals);
  }
  if (particle_density_) {
    particle_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
    particle_density_.print_result(xvals);
  }
  if (spin_corr_) {
    spin_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
    spin_corr_.print_result(xvals);
  }
  if (sc_corr_) {
    sc_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
    sc_corr_.print_result(xvals);
  }
  if (particle_density_) {
    particle_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
    particle_density_.print_result(xvals);
  }
  if (doublon_density_) {
    doublon_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
    doublon_density_.print_result(xvals);
  }
}

void ObservableSet::print_results(const double& xval) 
{
  if (num_xvars_ != 1) 
    throw std::invalid_argument("ObservableSet::print_result: 'x-vars' size mismatch");
  std::vector<double> xvals{xval};
  if (energy_) {
    energy_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_.print_result(xvals);
  }
  if (energy_grad_) {
    energy_grad_.print_heading(headstream_.rdbuf()->str(),xvars_);
    energy_grad_.print_result(xvals);
  }
  if (spin_corr_) {
    spin_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
    spin_corr_.print_result(xvals);
  }
  if (sc_corr_) {
    sc_corr_.print_heading(headstream_.rdbuf()->str(),xvars_);
    sc_corr_.print_result(xvals);
  }
  if (particle_density_) {
    particle_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
    particle_density_.print_result(xvals);
  }
  if (doublon_density_) {
    doublon_density_.print_heading(headstream_.rdbuf()->str(),xvars_);
    doublon_density_.print_result(xvals);
  }
}

void ObservableSet::MPI_send_results(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& proc, const int& msg_tag)
{
  if (energy_) energy_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (energy_grad_) energy_grad_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (spin_corr_) spin_corr_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (sc_corr_) sc_corr_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (sr_matrix_) sr_matrix_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (particle_density_) particle_density_.MPI_send_data(mpi_comm, proc, msg_tag);
  if (doublon_density_) doublon_density_.MPI_send_data(mpi_comm, proc, msg_tag);
}

void ObservableSet::MPI_recv_results(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& proc, const int& msg_tag)
{
  if (energy_) energy_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (energy_grad_) energy_grad_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (spin_corr_) spin_corr_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (sc_corr_) sc_corr_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (sr_matrix_) sr_matrix_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (particle_density_) particle_density_.MPI_add_data(mpi_comm, proc, msg_tag);
  if (doublon_density_) doublon_density_.MPI_add_data(mpi_comm, proc, msg_tag);
}

} // end namespace vmc


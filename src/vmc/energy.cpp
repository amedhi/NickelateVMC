/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:47:12
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 00:13:28
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./energy.h"

namespace vmc {

//-------------------ENERGY--------------------------------
void Energy::setup(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const bool& only_total)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  std::vector<std::string> elem_names;
  std::vector<std::string> term_names;
  //model.get_term_names(elem_names);
  //this->resize(elem_names.size(), elem_names);
  //this->set_have_total();

  // include the 'Total' as seperate component
  elem_names.push_back("Total");
  model.get_term_names(term_names);
  for (const auto& s : term_names) elem_names.push_back(s);
  config_value_.resize(elem_names.size());

  // if only total component is required
  only_total_ = only_total;
  if (only_total_) {
    elem_names.clear();
    elem_names.push_back("Energy");
    this->resize(elem_names.size(), elem_names);
  }
  else {
    this->resize(elem_names.size(), elem_names);
  }

  setup_done_ = true;
}
  
void Energy::measure(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const SysConfig& config, 
  const SiteDisorder& site_disorder)
{
  using op_id = model::op_id;
  // In the energy components below, 0-th component is for the total
  int component = 1;
  config_value_.setZero();
  // bond energies
  if (model.have_bondterm()) {
    Matrix matrix_elem(model.num_bondterms(),lattice.num_bond_types());
    matrix_elem.setZero();
    for (const auto& b : lattice.bonds()) {
      unsigned type = b.type();
      unsigned site_i = b.src_id();
      unsigned site_j = b.tgt_id();
      // matrix elements each term & bond type
      int bterm = 0;
      for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
        // compute only if non-zero coupling constant
        if (std::abs(it->coupling(type))>1.0E-8) {
          matrix_elem(bterm,type) += config.apply(it->qn_operator(),site_i,site_j,
            b.sign(),b.phase());
        }
        bterm++;
      }
    }
    // sum the contributions
    int bterm = 0;
    for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
      for (unsigned btype=0; btype<lattice.num_bond_types(); ++btype) {
        config_value_(component) += std::real(it->coupling(btype)*matrix_elem(bterm,btype));
        //std::cout << matrix_elem(bterm,btype) << "\n";
        //std::cout << btype << "  " << it->coupling(btype) << "\n"; getchar();
      }
      bterm++;
      component++;
    }
  }

  // site energies
  if (model.have_siteterm()) {
    Matrix matrix_elem(model.num_siteterms(),lattice.num_site_types());
    matrix_elem.setZero();
    Eigen::VectorXi hubbard_nd(lattice.num_site_types());
    hubbard_nd.setZero();
    // do it little differently
    int sterm = 0;
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (const auto& s : lattice.sites()) {
          hubbard_nd(s.type()) += config.apply_ni_dblon(s.id());
        }
      }
      else {
        for (const auto& s : lattice.sites()) {
          matrix_elem(sterm,s.type()) += config.apply(it->qn_operator(), s.id());
          //std::cout << config.apply(it->qn_operator(),site) << "\n"; getchar();
        }
      }
      sterm++;
    }
    sterm = 0;
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (int stype=0; stype<lattice.num_site_types(); ++stype) {
          config_value_(component) += std::real(it->coupling(stype))*hubbard_nd(stype);
        }
      }
      else {
        for (unsigned stype=0; stype<lattice.num_site_types(); ++stype) {
          config_value_(component) += std::real(it->coupling(stype)*matrix_elem(sterm,stype));
        }
      }
      sterm++;
      component++;
    }
  }

  // disorder term
  if (site_disorder) {
    //std::cout << "\ndisorder energy\n"; 
    double disorder_en = 0.0;
    for (const auto& s : lattice.sites()) {
      int n_i = config.apply(model::op::ni_sigma(), s.id());
      disorder_en += std::real(n_i * site_disorder.potential(s.id()));
      //std::cout <<"site= "<<site<<" ni= "<<n_i;
      //std::cout <<" V= "<<site_disorder.potential(site)<<"\n";
      //std::cout << "E+ = " << disorder_en << "\n"; 
    }
    config_value_(component) = disorder_en;
    //std::cout << "\ndisorder_en = " << disorder_en << "\n"; getchar();
  }

  // energy per site
  double total = config_value_.sum();

  // add to databin
  if (only_total_) {
    *this << total/num_sites_;
  }
  else {
    config_value_(0) = total;
    config_value_ /= num_sites_;
    *this << config_value_;
  }
  //std::cout << "config_energy = " << config_value_.transpose() << "\n"; 
  //getchar();
}

/*
void Energy::setup(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  std::vector<std::string> elem_names;
  std::vector<std::string> term_names;
  //model.get_term_names(elem_names);
  //this->resize(elem_names.size(), elem_names);
  //this->set_have_total();

  // include the 'Total' as seperate component
  model.get_term_names(term_names);
  elem_names.push_back("Total");
  for (const auto& s : term_names) elem_names.push_back(s);
  this->resize(elem_names.size(), elem_names);
  config_value_.resize(elem_names.size());
  setup_done_ = true;
}
  
void Energy::measure(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const SysConfig& config, 
  const SiteDisorder& site_disorder)
{
  using op_id = model::op_id;
  // In the energy components below, 0-th component is for the total
  int component = 1;
  config_value_.setZero();
  // bond energies
  if (model.have_bondterm()) {
    Matrix matrix_elem(model.num_bondterms(),lattice.num_bond_types());
    matrix_elem.setZero();
    for (const auto& b : lattice.bonds()) {
      unsigned type = b.type();
      unsigned site_i = b.src_id();
      unsigned site_j = b.tgt_id();
      // matrix elements each term & bond type
      int bterm = 0;
      for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
        matrix_elem(bterm,type) += config.apply(it->qn_operator(),site_i,site_j,
          b.sign(),b.phase());
        bterm++;
      }
    }
    // sum the contributions
    int bterm = 0;
    for (auto it=model.bondterms_begin(); it!=model.bondterms_end(); ++it) {
      for (unsigned btype=0; btype<lattice.num_bond_types(); ++btype) {
        config_value_(component) += std::real(it->coupling(btype)*matrix_elem(bterm,btype));
        //std::cout << matrix_elem(bterm,btype) << "\n";
        //std::cout << btype << "  " << it->coupling(btype) << "\n"; getchar();
      }
      bterm++;
      component++;
    }
  }

  // site energies
  if (model.have_siteterm()) {
    Matrix matrix_elem(model.num_siteterms(),lattice.num_site_types());
    matrix_elem.setZero();
    Eigen::VectorXi hubbard_nd(lattice.num_site_types());
    hubbard_nd.setZero();
    // do it little differently
    int sterm = 0;
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (const auto& s : lattice.sites()) {
          hubbard_nd(s.type()) += config.apply_ni_dblon(s.id());
        }
      }
      else {
        for (const auto& s : lattice.sites()) {
          matrix_elem(sterm,s.type()) += config.apply(it->qn_operator(), s.id());
          //std::cout << config.apply(it->qn_operator(),site) << "\n"; getchar();
        }
      }
      sterm++;
    }
    sterm = 0;
    for (auto it=model.siteterms_begin(); it!=model.siteterms_end(); ++it) {
      // special treatment for hubbard
      if (it->qn_operator().id()==op_id::niup_nidn) {
        for (int stype=0; stype<lattice.num_site_types(); ++stype) {
          config_value_(component) += std::real(it->coupling(stype))*hubbard_nd(stype);
        }
      }
      else {
        for (unsigned stype=0; stype<lattice.num_site_types(); ++stype) {
          config_value_(component) += std::real(it->coupling(stype)*matrix_elem(sterm,stype));
        }
      }
      sterm++;
      component++;
    }
  }

  // disorder term
  if (site_disorder) {
    //std::cout << "\ndisorder energy\n"; 
    double disorder_en = 0.0;
    for (const auto& s : lattice.sites()) {
      int n_i = config.apply(model::op::ni_sigma(), s.id());
      disorder_en += std::real(n_i * site_disorder.potential(s.id()));
      //std::cout <<"site= "<<site<<" ni= "<<n_i;
      //std::cout <<" V= "<<site_disorder.potential(site)<<"\n";
      //std::cout << "E+ = " << disorder_en << "\n"; 
    }
    config_value_(component) = disorder_en;
    //std::cout << "\ndisorder_en = " << disorder_en << "\n"; getchar();
  }

  // energy per site
  double total = config_value_.sum();
  config_value_(0) = total;
  config_value_ /= num_sites_;
  //std::cout << "config_energy = " << config_value_.transpose() << "\n"; 
  //getchar();

  // add to databin
  *this << config_value_;
}
*/

//-------------------ENERGY GRADIENT--------------------------------
void EnergyGradient::setup(const SysConfig& config, const int& sample_size)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_varp_ = config.num_varparms();
  this->resize(config.num_varparms(),config.varp_names());
  grad_logpsi_.resize(config.num_varparms());
  config_value_.resize(2*config.num_varparms());

  // for calculating gradients from batch data
  assert(sample_size>=0);
  batch_size_ = std::min(sample_size,500);
  bool no_error_bar = true;
  batch_gsum_.init("batch_gsum",2*config.num_varparms(),no_error_bar);
  batch_esum_.init("batch_esum",1,no_error_bar);

  setup_done_ = true;
}

void EnergyGradient::reset(void) 
{
  MC_Observable::reset(); 
  batch_gsum_.reset(); 
  batch_esum_.reset();
}

void EnergyGradient::reset_batch_limit(const int& sample_size)
{
  assert(sample_size>=0);
  batch_size_ = std::min(sample_size,500);
}
  
void EnergyGradient::measure(const SysConfig& config, const double& config_energy)
{
  config.get_grad_logpsi(grad_logpsi_);
  int n = 0;
  for (int i=0; i<num_varp_; ++i) {
    config_value_[n] = config_energy*grad_logpsi_[i];
    config_value_[n+1] = grad_logpsi_[i];
    n += 2;
  }
  batch_gsum_ << config_value_;
  batch_esum_ << config_energy;

  // compute gradient from the batch sum 
  if (batch_esum_.num_samples()>=batch_size_) {
    mcdata::data_t grad(num_varp_);
    double Emean = batch_esum_.mean();
    config_value_ = batch_gsum_.mean_data(); 
    int n = 0;
    for (int i=0; i<num_varp_; ++i) {
      grad(i) = 2.0*(config_value_[n] - Emean*config_value_[n+1]);
      n += 2;
    }
    // add sample
    *this << grad;
    // reset partial sums
    batch_gsum_.reset();
    batch_esum_.reset();
  }
}

void EnergyGradient::finalize(void)
{
  // nothing to do
}

/*
void EnergyGradient::measure(const SysConfig& config, const double& config_energy)
{
  config.get_grad_logpsi(grad_logpsi_);
  unsigned n = 0;
  for (unsigned i=0; i<num_varp_; ++i) {
    config_value_[n] = config_energy*grad_logpsi_[i];
    config_value_[n+1] = grad_logpsi_[i];
    n += 2;
  }
  grad_terms_ << config_value_;
  total_en_ << config_energy;
}

void EnergyGradient::finalize(void)
{
  double mean_energy = total_en_.mean();
  config_value_ = grad_terms_.mean_data(); 
  mcdata::data_t energy_grad(num_varp_);
  int n = 0;
  for (int i=0; i<num_varp_; ++i) {
    energy_grad(i) = 2.0*(config_value_[n] - mean_energy*config_value_[n+1]);
    n += 2;
  }
  //-------------------------------------------------
  //std::cout << mean_energy << "\n";
  //std::cout << config_value_.transpose() << "\n";
  //for (int i=0; i<energy_grad_.size(); ++i) std::cout << energy_grad_[i] << "\n";
  //-------------------------------------------------
  *this << energy_grad;
  // can we reset these now? not reseting...
  //grad_terms_.reset();
  //total_en_.reset();
}
*/

/*
#ifdef HAVE_BOOST_MPI
void EnergyGradient::MPI_send_data(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& p, const int& msg_tag)
{
  mcdata::MC_Observable::MPI_send_data(mpi_comm, p, msg_tag);
  grad_terms_.MPI_send_data(mpi_comm, p, msg_tag);
  total_en_.MPI_send_data(mpi_comm, p, msg_tag);
}

void EnergyGradient::MPI_add_data(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& p, const int& msg_tag)
{
  mcdata::MC_Observable::MPI_add_data(mpi_comm, p, msg_tag);
  grad_terms_.MPI_add_data(mpi_comm, p, msg_tag);
  total_en_.MPI_add_data(mpi_comm, p, msg_tag);
}
#else
void EnergyGradient::MPI_send_data(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& p, const int& msg_tag)
{
  mcdata::MC_Observable::MPI_send_data(mpi_comm, p, msg_tag);
}
void EnergyGradient::MPI_add_data(const mpi::mpi_communicator& mpi_comm, 
  const mpi::proc& p, const int& msg_tag)
{
  mcdata::MC_Observable::MPI_add_data(mpi_comm, p, msg_tag);
}
#endif
*/


//-------------------SR_Matrix (Stochastic Reconfiguration)---------------
void SR_Matrix::setup(const lattice::Lattice& lattice, const SysConfig& config,
  const int& sample_size)
{
  MC_Observable::switch_on();
  if (setup_done_) return;
  num_sites_ = lattice.num_sites();
  num_varp_ = config.num_varparms();
  // '\del(ln(psi))' plus upper triangular part of the sr_matrix 
  unsigned num_elems = num_varp_ + num_varp_*(num_varp_+1)/2;
  this->resize(num_elems);
  this->error_bar_off();
  config_value_.resize(num_elems);
  upper_elems_.resize(num_elems);

  // for calculating gradients from batch data
  assert(sample_size>=0);
  batch_size_ = std::min(sample_size,500);
  batch_sum_.resize(num_elems);
  batch_sum_.error_bar_off();
  setup_done_ = true;
}

void SR_Matrix::reset(void) 
{
  MC_Observable::reset(); 
  batch_sum_.reset(); 
}

void SR_Matrix::reset_batch_limit(const int& sample_size)
{
  assert(sample_size>=0);
  batch_size_ = std::min(sample_size,500);
}
  

void SR_Matrix::measure(const RealVector& grad_logpsi)
{
  assert(grad_logpsi.size()==num_varp_);
 /*--------------------------------------------------
  * Construct the upper triangular part of Smatrix, 
  * and flatten it into a vector
  *--------------------------------------------------*/
  // operator 'del(ln(psi))' terms
  for (unsigned i=0; i<num_varp_; ++i) config_value_[i] = grad_logpsi[i];

  unsigned k = num_varp_;
  for (unsigned i=0; i<num_varp_; ++i) {
    double x = grad_logpsi[i];
    for (unsigned j=i; j<num_varp_; ++j) {
      double y = grad_logpsi[j];
      config_value_[k] = x * y;
      ++k;
    }
  }
  // add to batch data
  batch_sum_ << config_value_;

 /*--------------------------------------------------
  * Compute the SR matrix (upper triangular) elements 
  * from the batch data 
  *--------------------------------------------------*/
  if (batch_sum_.num_samples()>=batch_size_) {
    config_value_ = batch_sum_.mean_data();
    unsigned k = num_varp_;
    int n = 0;
    for (unsigned i=0; i<num_varp_; ++i) {
      double x = config_value_[i];
      for (unsigned j=i; j<num_varp_; ++j) {
        double y = config_value_[j];
        upper_elems_[n] = (config_value_[k] - x*y)/num_sites_;
        k++;
        n++;
      }
    }
    // add sample
    *this << upper_elems_;
    // reset partial sums
    batch_sum_.reset();
  }
}

void SR_Matrix::get_matrix(Eigen::MatrixXd& sr_matrix) const
{
  // construct the matrix from the mean elemenet values
  upper_elems_ = MC_Observable::mean_data();
  int n = 0;
  for (unsigned i=0; i<num_varp_; ++i) {
    for (unsigned j=i; j<num_varp_; ++j) {
      sr_matrix(i,j) = upper_elems_[n];
      sr_matrix(j,i) = sr_matrix(i,j);
      n++;
    }
  }
}

/*
void SR_Matrix::measure(const RealVector& grad_logpsi)
{
  assert(grad_logpsi.size()==num_varp_);
  // operator 'del(ln(psi))' terms
  for (unsigned i=0; i<num_varp_; ++i) config_value_[i] = grad_logpsi[i];
  // flatten the upper triangular part to a vector
  unsigned k = num_varp_;
  for (unsigned i=0; i<num_varp_; ++i) {
    double x = grad_logpsi[i];
    for (unsigned j=i; j<num_varp_; ++j) {
      double y = grad_logpsi[j];
      config_value_[k] = x * y;
      ++k;
    }
  }
  *this << config_value_;
}

void SR_Matrix::get_matrix(Eigen::MatrixXd& sr_matrix) const
{
  // 'config_value' is used as temporary storage
  config_value_ = MC_Observable::mean_data();
  unsigned k = num_varp_;
  for (unsigned i=0; i<num_varp_; ++i) {
    double x = config_value_[i];
    for (unsigned j=i; j<num_varp_; ++j) {
      double y = config_value_[j];
      sr_matrix(i,j) = (config_value_[k] - x*y)/num_sites_;
      sr_matrix(j,i) = sr_matrix(i,j);
      ++k;
    }
  }
}
*/


} // end namespace vmc

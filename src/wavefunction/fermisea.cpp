/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-02-20 12:21:42
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2023-06-20 17:51:17
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <numeric>
#include "./fermisea.h"
#include <fstream>
#include <iomanip>
#include <cassert>

namespace var {

Fermisea::Fermisea(const MF_Order::order_t& order, const input::Parameters& inputs, 
  const lattice::Lattice& lattice, const model::Hamiltonian& model) 
  : GroundState(order, MF_Order::pairing_t::null)
{
  init(inputs, lattice, model);
}

int Fermisea::init(const input::Parameters& inputs, const lattice::Lattice& lattice,
  const model::Hamiltonian& model)
{
  name_ = "Fermisea";
  lattice_id_ = lattice.id();
  // sites & bonds
  num_sites_ = lattice.num_sites();
  num_bonds_ = lattice.num_bonds();

  // particle number
  set_nonmagnetic(true);
  set_particle_num(inputs);

  // bloch basis
  blochbasis_.construct(lattice);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(lattice);

  // build MF Hamiltonian
  varparms_.clear();
  mf_model_.init(lattice);
  std::string name;
  double defval, lb, ub, dh;
  using namespace model;
  model::CouplingConstant cc;

  // SC form
  correlation_pairs().clear();
  const std::pair<int,int> anypair{-1,-1};
  using order_t = MF_Order::order_t;
  int info;

  // Ground state degeneracy warning
  degeneracy_warning_ = inputs.set_value("warn_degeneracy", true, info);

  /*-------------------------------------------------
  * Chemical Potential:
  * Unless 'mu' is to be taken as variational parameter,
  * it is taken to be '0', by default. 
  *--------------------------------------------------*/
  mu_variational_ = inputs.set_value("mu_variational", false, info);


  if (lattice.id()==lattice::lattice_id::SQUARE_NNN) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=1.0, inputs);
    cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
    if (mu_variational_) {
      mf_model_.add_parameter(name="mu",defval=0.0,inputs);
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      defval = mf_model_.get_parameter_value("mu");
      varparms_.add("mu",defval, lb=-5.0, ub=+5.0, dh=0.1);
    }
  }

  else if (lattice.id()==lattice::lattice_id::SQUARE_2SITE) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_af", defval=1.0, inputs);

    // bond operator terms
    cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

    // site operator terms
    cc.create(2);
    cc.add_type(0, "-delta_af");
    cc.add_type(1, "delta_af");
    mf_model_.add_siteterm(name="ni_up", cc, op::ni_up());
    cc.create(2);
    cc.add_type(0, "delta_af");
    cc.add_type(1, "-delta_af");
    mf_model_.add_siteterm(name="ni_dn", cc, op::ni_dn());

    // variational parameters
    defval = mf_model_.get_parameter_value("delta_af");
    varparms_.add("delta_af",defval,lb=0,ub=2.0,dh=0.02);
    if (mu_variational_) {
      mf_model_.add_parameter(name="mu",defval=0.0,inputs);
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      defval = mf_model_.get_parameter_value("mu");
      varparms_.add("mu", defval, lb=-5.0, ub=+5.0, dh=0.1);
    }
  }

  else if (lattice.id()==lattice::lattice_id::SQUARE_4SITE) {
    if (model.id()==model::model_id::HUBBARD_IONIC || 
        model.id()== model::model_id::TJ_IONIC) {
      bool mu_default = inputs.set_value("mu_default", false);
      mf_model_.add_parameter(name="tv", defval=1.0, inputs);
      mf_model_.add_parameter(name="tpv", defval=1.0, inputs);
      mf_model_.add_parameter(name="dAF", defval=0.0, inputs);
      mf_model_.add_parameter(name="muA", defval=0.0, inputs);
      mf_model_.add_parameter(name="muB", defval=0.0, inputs);

      // A-B hopping term
      cc.create(6);
      cc.add_type(0,"-tv");
      cc.add_type(1,"-tv");
      cc.add_type(2,"0");
      cc.add_type(3,"0");
      cc.add_type(4,"0");
      cc.add_type(5,"0");
      mf_model_.add_bondterm(name="hop-AB", cc, op::spin_hop());

      // A-A hopping term
      cc.create(6);
      cc.add_type(0,"0");
      cc.add_type(1,"0");
      cc.add_type(2,"-tpv");
      cc.add_type(3,"-tpv");
      cc.add_type(4,"0");
      cc.add_type(5,"0");
      mf_model_.add_bondterm(name="hop-AA", cc, op::spin_hop());

      // B-B hopping term
      cc.create(6);
      cc.add_type(0,"0");
      cc.add_type(1,"0");
      cc.add_type(2,"0");
      cc.add_type(3,"0");
      cc.add_type(4,"-tpv");
      cc.add_type(5,"-tpv");
      mf_model_.add_bondterm(name="hop-BB", cc, op::spin_hop());

      // Site terms: AF-order & chemical potential
      cc.create(2);
      cc.add_type(0, "-(muA+dAF)");
      cc.add_type(1, "-(muB-dAF)");
      mf_model_.add_siteterm(name="ni_up", cc, op::ni_up());
      cc.create(2);
      cc.add_type(0, "-(muA-dAF)");
      cc.add_type(1, "-(muB+dAF)");
      mf_model_.add_siteterm(name="ni_dn", cc, op::ni_dn());

      correlation_pairs().push_back({0,0});
      correlation_pairs().push_back({0,1});
      correlation_pairs().push_back({2,2});
      correlation_pairs().push_back({2,3});
      correlation_pairs().push_back({4,4});
      correlation_pairs().push_back({4,5});

      // variational parameters
      defval = mf_model_.get_parameter_value("dAF");
      varparms_.add("dAF", defval,lb=0.0,ub=+5.0,dh=0.1);

      // chemical potential
      mu_variational_ = mu_default; // for this model
      if (mu_variational_) {
        defval = mf_model_.get_parameter_value("muA");
        varparms_.add("muA", defval,lb=-10.0,ub=+10.0,dh=0.1);
        defval = mf_model_.get_parameter_value("muB");
        varparms_.add("muB", defval,lb=-10.0,ub=+10.0,dh=0.1);
      }

      // hopping integrals as variational parameters
      if (inputs.set_value("tv_variational", true)) {
        defval = mf_model_.get_parameter_value("tv");
        varparms_.add("tv", defval,lb=0.1,ub=5.0,dh=0.02);
        defval = mf_model_.get_parameter_value("tpv");
        varparms_.add("tpv", defval,lb=0.0,ub=2.0,dh=0.02);
      }
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for the model");
    }
  }

  else if (lattice.id()==lattice::lattice_id::CHAIN_2SITE) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);

    // bond operator terms
    cc = CouplingConstant({0,"-t"},{1,"-t"});
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

    if (order()==order_t::AF) {
      order_name_ = "AF";
      mf_model_.add_parameter(name="delta_af", defval=1.0, inputs);
      // site operator terms
      cc.create(2);
      cc.add_type(0, "-delta_af");
      cc.add_type(1, "delta_af");
      mf_model_.add_siteterm(name="ni_up", cc, op::ni_up());
      cc.create(2);
      cc.add_type(0, "delta_af");
      cc.add_type(1, "-delta_af");
      mf_model_.add_siteterm(name="ni_dn", cc, op::ni_dn());

      // variational parameters
      defval = mf_model_.get_parameter_value("delta_af");
      varparms_.add("delta_af",defval,lb=0,ub=2.0,dh=0.02);
      if (mu_variational_) {
        mf_model_.add_parameter(name="mu",defval=0.0,inputs);
        mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
        defval = mf_model_.get_parameter_value("mu");
        varparms_.add("mu", defval, lb=-5.0, ub=+5.0, dh=0.1);
      }
    }
  }

  else if (lattice.id()==lattice::lattice_id::SQUARE_STRIPE) {
    mf_model_.add_parameter(name="tx", defval=1.0, inputs);
    mf_model_.add_parameter(name="ty1", defval=1.0, inputs);
    mf_model_.add_parameter(name="ty2", defval=1.0, inputs);
    mf_model_.add_parameter(name="lambda", defval=8, inputs);
    // variational parameter
    mf_model_.add_parameter(name="dCDW", defval=0.0, inputs);
    mf_model_.add_parameter(name="dSDW", defval=0.0, inputs);

    // hopping term 
    cc.create(48);
    for (int i=0; i<16; ++i) {
      cc.add_type(i, "-tx");
    }
    for (int i=16; i<32; ++i) {
      cc.add_type(i, "-ty1");
    }
    for (int i=32; i<48; ++i) {
      cc.add_type(i, "-ty2");
    }
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

    std::ostringstream strval;
    strval.precision(12);
    strval.setf(std::ios_base::fixed);
    int lambda = mf_model_.get_parameter_value("lambda");
    // check commensurate or not
    double Q;
    if (lambda > 0) {
      int m = 16/lambda;
      if (16 != m*lambda) {
        std::cout<<">> alert: BCS_State::init: incommensurate `lambda`\n";
      }
      Q = two_pi()/lambda;
    }
    else {
      Q = 0.0;
    }

    // CDW term
    cc.create(32);
    for (int n=0; n<32; ++n) {
      double ampl;
      if (lambda > 0) ampl = std::cos(Q*(n+1-0.5));
      else ampl = 1.0;
      strval << ampl << "*dCDW";
      cc.add_type(n, strval.str());
      //std::cout << "ampl["<<n<<"] = "<< strval.str() << "\n";
      strval.str(std::string()); // clear
    }
    mf_model_.add_siteterm(name="CDW", cc, op::ni_sigma());

    // SDW term (UP spin)
    cc.create(32);
    strval.str(std::string()); // clear
    int sign = 1;
    for (int n=0; n<32; ++n) {
      if (n == 16) sign = -1;
      double ampl;
      if (lambda > 0) ampl = sign*std::sin(0.5*Q*(n+1-0.5));
      else ampl = sign;
      sign *= -1;
      strval << ampl << "*dSDW";
      cc.add_type(n, strval.str());
      //std::cout << "ampl["<<n<<"] = "<< strval.str() << "\n";
      strval.str(std::string()); // clear
    }
    mf_model_.add_siteterm(name="SDW_UP", cc, op::ni_up());

    // SDW term (DN spin)
    cc.create(32);
    strval.str(std::string()); // clear
    sign = -1;
    for (int n=0; n<32; ++n) {
      if (n == 16) sign = 1;
      double ampl;
      if (lambda > 0) ampl = sign*std::sin(0.5*Q*(n+1-0.5));
      else ampl = sign;
      sign *= -1;
      strval << ampl << "*dSDW";
      cc.add_type(n, strval.str());
      //std::cout << "ampl["<<n<<"] = "<< strval.str() << "\n";
      strval.str(std::string()); // clear
    }
    mf_model_.add_siteterm(name="SDW_DN", cc, op::ni_dn());

    // variational parameters
    defval = mf_model_.get_parameter_value("dCDW");
    varparms_.add("dCDW", defval,lb=0,ub=2.0,dh=0.01);
    defval = mf_model_.get_parameter_value("dSDW");
    varparms_.add("dSDW", defval,lb=0,ub=2.0,dh=0.01);
  }

  else if (lattice.id()==lattice::lattice_id::NICKELATE_2L) {
    mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=0.0, inputs);
    mf_model_.add_parameter(name="th", defval=0.0, inputs);

    // bond operators
    cc.create(7);
    cc.add_type(0, "-t");
    cc.add_type(1, "-t");
    cc.add_type(2, "-tp");
    cc.add_type(3, "-t");
    cc.add_type(4, "-t");
    cc.add_type(5, "-tp");
    cc.add_type(6, "-th");
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
    // site operators
    cc.create(2);
    cc.add_type(0, "0");
    cc.add_type(1, "e_R");
    mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());
    if (mu_variational_) {
      mf_model_.add_parameter(name="mu",defval=0.0,inputs);
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      defval = mf_model_.get_parameter_value("mu");
      varparms_.add("mu", defval, lb=-5.0, ub=+5.0, dh=0.1);
    }
  }

  else {
    throw std::range_error("Fermisea::Fermisea: state undefined for the lattice");
  }

  // finalize MF Hamiltonian
  num_varparms_ = varparms_.size();
  mf_model_.finalize(lattice);

  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  work2_.resize(num_sites_,num_sites_);
  phi_k_.resize(num_kpoints_);
  work_k_.resize(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
    work_k_[k].resize(kblock_dim_,kblock_dim_);
  } 

  // for calculating SC correlations 
  rmax_ =  lattice.size1()/2+1;
  alpha_ = lattice.basis_vector_a1();
  beta_ = lattice.basis_vector_a2();
  R_list_.resize(rmax_);
  corr_aa_.resize(kblock_dim_,rmax_);
  corr_ab_.resize(kblock_dim_,rmax_);
  corr_fs_.resize(kblock_dim_,rmax_);
  for (int r=0; r<rmax_; ++r) {
    R_list_[r] = r*alpha_;
    //std::cout << r << "  " << R_list_[r].transpose() << "\n";
  }
  //getchar();

  return 0;
}


void Fermisea::update(const lattice::Lattice& lattice)
{
  // update for change in lattice BC (same structure & size)
  // bloch basis
  blochbasis_.construct(lattice);
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(lattice);
}

std::string Fermisea::info_str(void) const
{
  std::ostringstream info;
  //info << "# Ground State: '"<<name_<<"'\n";
  info << "# Ground State: '"<<name_<<" ("<<order_name_<<")'\n";
  info << "# Hole doping = "<<hole_doping()<<"\n";
  info << "# Particles = "<< num_upspins()+num_dnspins();
  info << " (Nup = "<<num_upspins()<<", Ndn="<<num_dnspins()<<")\n";
  info.precision(6);
  return info.str();
}

void Fermisea::update(const input::Parameters& inputs)
{
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
  // update MF model
  mf_model_.update(inputs);
  // update variational parameters
  for (auto& p : varparms_) 
    p.change_value(mf_model_.get_parameter_value(p.name()));
}

void Fermisea::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
}

void Fermisea::get_wf_amplitudes(Matrix& psi) 
{
  construct_groundstate();
  get_pair_amplitudes_sitebasis2(psi);
  //get_pair_amplitudes(phi_k_);
  //get_pair_amplitudes_sitebasis(phi_k_, psi);
  /*std::cout << "Fermisea::update: SC correlations\n";
  get_sc_correlation();
  getchar();
  */
}

void Fermisea::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  // numerical gradient
  int i=0; 
  for (const auto& p : varparms_) {
    double h = p.diff_h();
    double x = p.value();
    double x_fwd = x+h; 
    if (x_fwd > p.ubound()) x_fwd = p.ubound();
    mf_model_.update_parameter(p.name(), x_fwd);
    mf_model_.update_terms();
    construct_groundstate();
    get_pair_amplitudes_sitebasis2(psi_gradient[i]);
    double x_bwd = x-h; 
    if (x_bwd < p.lbound()) x_bwd = p.lbound();
    mf_model_.update_parameter(p.name(), x_bwd);
    mf_model_.update_terms();
    construct_groundstate();
    get_pair_amplitudes_sitebasis2(work2_);
    // restore model to original state
    mf_model_.update_parameter(p.name(), x);
    mf_model_.update_terms();
    construct_groundstate();
    // finite difference gradients
    double inv_2h = 1.0/(x_fwd-x_bwd);
    psi_gradient[i] -= work2_;
    psi_gradient[i] *= inv_2h;
    ++i;
  }
}

/*
void Fermisea::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  // numerical gradient
  int i=0; 
  for (const auto& p : varparms_) {
    double h = p.diff_h();
    double x = p.value();
    double x_fwd = x+h; 
    if (x_fwd > p.ubound()) x_fwd = p.ubound();
    mf_model_.update_parameter(p.name(), x_fwd);
    mf_model_.update_terms();
    get_pair_amplitudes(phi_k_);
    double x_bwd = x-h; 
    if (x_bwd < p.lbound()) x_bwd = p.lbound();
    mf_model_.update_parameter(p.name(), x_bwd);
    mf_model_.update_terms();
    get_pair_amplitudes(phi_k_);
    // model to original state
    mf_model_.update_parameter(p.name(), x);
    mf_model_.update_terms();
    //double inv_2h = 0.5/h;
    double inv_2h = 1.0/(x_fwd-x_bwd);
    for (int k=0; k<num_kpoints_; ++k) {
      phi_k_[k] -= work_k_[k];
      phi_k_[k] *= inv_2h;
    }
    //std::cout << phi_k_[0] << "\n"; getchar();
    // phi_k_ is now the derivative wrt i-th parameter
    // wave function gradients
    get_pair_amplitudes_sitebasis(phi_k_, psi_gradient[i]);
    //get_pair_amplitudes_sitebasis(psi_gradient[i]);
    //std::cout << psi_gradient[i] << "\n"; getchar();
    ++i;
  }
}
*/

void Fermisea::get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k) 
{
  for (int i=0; i<kshells_up_.size(); ++i) {
    int k = kshells_up_[i].k;
    int m = kshells_up_[i].nmax+1;
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_dn.compute(mf_model_.quadratic_spindn_block());
    phi_k[k] = es_k_up.eigenvectors().block(0,0,kblock_dim_,m)
      		 * es_minusk_dn.eigenvectors().transpose().block(0,0,m,kblock_dim_);
    //std::cout << "kvec = "<< kvec.transpose() << "\n"; 
    //std::cout << mf_model_.quadratic_spinup_block() << "\n"; 
    //std::cout << phi_k[k] << "\n"; getchar();
  }
}

void Fermisea::get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
  // psi = FTU_ * PHI_K * conjugate(transpose(FTU_))
  // PHI_K is block diagonal (k-th block is phi_k) 
  int p = 0;
  for (int i=0; i<num_kpoints_; ++i) {
    int q = 0;
    for (int j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (int ks=0; ks<kshells_up_.size(); ++ks) {
        int k = kshells_up_[ks].k;
        work_ += FTU_(i,k) * phi_k[k] * std::conj(FTU_(j,k));
      }
      // std::cout << work_ << "\n"; getchar();
      // copy transformed block
      //psi.block(p,q,kblock_dim_,kblock_dim_) = 
      for (int m=0; m<kblock_dim_; ++m) {
        for (int n=0; n<kblock_dim_; ++n) {
          psi(p+m,q+n) = ampl_part(work_(m,n));
        }
      }
      q += kblock_dim_;
    }
    p += kblock_dim_;
  }
  //std::cout << psi << "\n";
  //getchar();
  /*
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      std::cout << "psi["<<i<<","<<j<<"] = "<<psi(i,j)<<"\n";
      getchar();
    }
  }*/

}

void Fermisea::get_pair_amplitudes_sitebasis2(Matrix& psi)
{
  ComplexMatrix psi_up(num_sites(),num_upspins());
  int n = 0;
  for (int p=0; p<kshells_up_.size(); ++p) {
    int k = kshells_up_[p].k;
    int nband = kshells_up_[p].nmax+1;
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_spinup_block());

    int m = 0;
    for (int I=0; I<num_kpoints_; ++I) {
      psi_up.block(m,n,kblock_dim_,nband) 
        = es_k_up.eigenvectors().block(0,0,kblock_dim_,nband)*FTU_(I,k);
      m += kblock_dim_;
    }
    n += nband;
  }

  ComplexMatrix psi_dn(num_sites(), num_dnspins());
  n = 0;
  for (int p=0; p<kshells_dn_.size(); ++p) {
    int k = kshells_dn_[p].k;
    int nband = kshells_dn_[p].nmax+1;
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_minusk_dn.compute(mf_model_.quadratic_spindn_block());
    int m = 0;
    for (int I=0; I<num_kpoints_; ++I) {
      psi_dn.block(m,n,kblock_dim_,nband) 
        = es_minusk_dn.eigenvectors().block(0,0,kblock_dim_,nband)*FTU_(I,k);
      m += kblock_dim_;
    }
    n += nband;
  }

  // final matrix
#ifdef REAL_WAVEFUNCTION
  ComplexMatrix psi_cmplx = psi_up * psi_dn.transpose();
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      psi(i,j) = ampl_part(psi_cmplx(i,j));
      //std::cout << "psi["<<i<<","<<j<<"] = "<<psi(i,j)<<"\n";
      //getchar();
    }
  }
#else
  psi = psi_up * psi_dn.transpose();
#endif

  //exit(0);

  //std::cout << psi << "\n";
  //getchar();
  /*
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      std::cout << "psi["<<i<<","<<j<<"] = "<<psi(i,j)<<"\n";
      getchar();
    }
  }*/

}

double Fermisea::get_mf_energy(void)
{
  double mf_energy = 0.0;
  return mf_energy/num_sites_;
}

void Fermisea::construct_groundstate(void)
{
  const int UP = 1;
  const int DN = -1;
  int num_upspin = num_upspins();
  int num_dnspin = num_dnspins();
  int num_particle = num_spins();
  assert(num_particle == (num_upspin+num_dnspin));

  // for checking
  std::ios state(NULL);
  state.copyfmt(std::cout);
  std::cout<<std::fixed<<std::setprecision(6)<<std::right;

  std::vector<std::tuple<int,int,int>> qn_list; // list of (k,n,sigma)
  std::vector<double> ek;
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);

    // spin-UP block
    es_k_up.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
    ek.insert(ek.end(),es_k_up.eigenvalues().data(),
      es_k_up.eigenvalues().data()+kblock_dim_);
    //std::cout << kvec.transpose() << " " << es_k_up.eigenvalues() << "\n"; getchar();
    for (int n=0; n<kblock_dim_; ++n) {
      qn_list.push_back({k, n, UP});
    }

    // spin-DN block
    es_minusk_dn.compute(mf_model_.quadratic_spindn_block(), Eigen::EigenvaluesOnly);
    ek.insert(ek.end(),es_minusk_dn.eigenvalues().data(),
      es_minusk_dn.eigenvalues().data()+kblock_dim_);
    //std::cout << kvec.transpose() << " " << es_k_up.eigenvalues() << "\n"; getchar();
    for (int n=0; n<kblock_dim_; ++n) {
      qn_list.push_back({k, n, DN});
    }
  }

  // Indices in the original ek-array is sorted according to increasing ek
  std::vector<int> idx(ek.size());
  std::iota(idx.begin(),idx.end(),0);
  std::sort(idx.begin(),idx.end(),[&ek](const int& i1, const int& i2) 
    {return ek[i1] < ek[i2];});
  /* 
  for (int i=0; i<ek.size(); ++i) {
    int k, n, s;
    std::tie(k,n,s) = qn_list[idx[i]];
    std::cout<<"ek[idx["<<std::setw(4)<<i<<"]], k  n  s = "<<std::setw(10)<<ek[idx[i]];
    std::cout<<std::setw(4)<<k<<std::setw(4)<<n<<std::setw(4)<<s<<"\n";
  }
  getchar();
  */
    
  // mean energy
  /*double e0 = 0.0;
  for (int i=0; i<num_particles; ++i) {
    e0 += ek[idx[i]];
  }
  e0 = e0 / num_sites_;
  std::cout << "e0 = " << e0 << "\n"; getchar();
  */

  // check for degeneracy 
  double degeneracy_tol = 1.0E-12;
  int top_filled_level = num_particle-1;
  int num_degen_state = 0;
  int num_valence_particle = 0;
  fermi_energy_ = ek[idx[top_filled_level]];
  total_energy_ = 0.0;
  for (int i=0; i<=top_filled_level; ++i) {
    total_energy_ += ek[idx[i]];
  }
  total_energy_ = total_energy_/num_sites_;
  //std::cout << "Total KE = " << total_energy_ << "\n";
  //std::cout << "fermi energy = " << fermi_energy_ << "\n";

  /*
  std::cout << "filled states\n";
  for (int i=0; i<num_particle; ++i) {
    int state = idx[i]; 
    std::tie(k,n,s) = qn_list[state];
    std::cout << "ek[idx["<<i<<"]] = "<<ek[idx[i]]<<" "<<k<<" "<<n<<" "<<s<<"\n";
  }
  std::cout << "Empty states\n";
  for (int i=num_particle; i<ek.size(); ++i) {
    int state = idx[i]; 
    std::tie(k,n,s) = qn_list[state];
    std::cout << "ek[idx["<<i<<"]] = "<<ek[idx[i]]<<" "<<k<<" "<<n<<" "<<s<<"\n";
  }
  getchar();
  */


  // look upward in energy
  for (int i=top_filled_level+1; i<ek.size(); ++i) {
    if (std::abs(fermi_energy_-ek[idx[i]])>degeneracy_tol) break;
    num_degen_state++;
  }
  // look downward in energy (correct to count from 'top_filled_level')
  if (num_degen_state>0) {
    for (int i=top_filled_level; i>=0; --i) {
      if (std::abs(fermi_energy_ - ek[idx[i]])>degeneracy_tol) break;
      num_degen_state++;
      num_valence_particle++;
    }
    // warn
    if (degeneracy_warning_) {
      std::cout << " >> warning! Groundstate degeneracy: " << num_valence_particle <<
        " particles in " << num_degen_state << " states" << "\n";
    }
  }

  /* 
    Filled k-shells. A k-shell is a group of energy levels having same 
    value of quantum number k.
  */

  // core particles
  int num_core_particle = num_particle - num_valence_particle;
  /*
  std::cout << "total_particles = " << num_particle << "\n";
  std::cout << "core_particles = " << num_core_particle << "\n";
  std::cout << "valence_particles = " << num_valence_particle << "\n";
  getchar();
  */

  // find 'nmax' values of filled k-shells
  std::vector<int> upshell_nmax(num_kpoints_);
  std::vector<int> dnshell_nmax(num_kpoints_);

  for (auto& elem : upshell_nmax) elem = -1; // invalid default value
  for (auto& elem : dnshell_nmax) elem = -1; // invalid default value
  int upspin_count = 0; 
  int dnspin_count = 0; 
  int k, n, s;
  for (int i=0; i<num_core_particle; ++i) {
    int state = idx[i]; 
    std::tie(k,n,s) = qn_list[state];
    if (s == UP) {
      upspin_count++;
      if (upshell_nmax[k] < n) upshell_nmax[k] = n;
    }
    if (s == DN) {
      dnspin_count++;
      if (dnshell_nmax[k] < n) dnshell_nmax[k] = n;
    }
  }

  if (upspin_count != dnspin_count) {
    std::cout << "upspin_count = " << upspin_count << "\n";
    std::cout << "dnspin_count = " << dnspin_count << "\n";
    throw std::logic_error("* Fermisea::construct_groundstate:: consistency check-1 failed!");
  }

  // order the degen states in order of increasing 'band' index, then 'k' index
  int s1 = num_core_particle;
  int s2 = num_core_particle+num_degen_state;
  for (int i=s1; i<s2; ++i) {
    int n1 = std::get<1>(qn_list[idx[i]]);
    int k1 = std::get<0>(qn_list[idx[i]]);
    double p1 = blochbasis_.kvector(k1).norm();
    for (int j=i+1; j<s2; ++j) {
      int n2 = std::get<1>(qn_list[idx[j]]);
      int k2 = std::get<0>(qn_list[idx[j]]);
      double p2 = blochbasis_.kvector(k2).norm();
      if ((n1>n2) || ((n1==n2) && (p1>p2))) {
        // swap
        int jdx = idx[j];
        idx[j] = idx[i];
        idx[i] = jdx;
      }
    }
  }
  /*
  std::cout << "Degenerate states\n";
  for (int i=s1; i<s2; ++i) {
    int k, n, s;
    std::tie(k,n,s) = qn_list[idx[i]];
    std::cout<<"ek[idx["<<i<<"]], k  n  s = "<<std::setw(10)<<ek[idx[i]];
    std::cout<<std::setw(4)<<k<<std::setw(4)<<n<<std::setw(4)<<s<<"\n";
  }
  getchar();
  */

  // fill up the valence particles
  for (int i=num_core_particle; i<num_particle; ++i) {
    int k, n, s; 
    std::tie(k,n,s) = qn_list[idx[i]];
    if (s==UP) {
      upspin_count++;
      if (upshell_nmax[k] < n) upshell_nmax[k] = n;
    }
    if (s==DN) {
      dnspin_count++;
      if (dnshell_nmax[k] < n) dnshell_nmax[k] = n;
    }
  }
  if ((upspin_count+dnspin_count) != num_particle) {
    throw std::logic_error("* Fermisea::construct_groundstate:: consistency check-2 failed!");
  }

  // currently, 'Sz /= 0' case is not implemented
  if (upspin_count != dnspin_count) {
    throw std::logic_error("* Fermisea::construct_groundstate:: Found Sz /= 0 state!");
  }

  // store the filled k-shells
  upspin_count = 0;
  kshells_up_.clear();
  for (int k=0; k<num_kpoints_; ++k) {
    int nmax = upshell_nmax[k];
    upspin_count += (nmax+1);
    if (nmax != -1) kshells_up_.push_back({k,0,nmax});
  }
  dnspin_count = 0;
  kshells_dn_.clear();
  for (int k=0; k<num_kpoints_; ++k) {
    int nmax = dnshell_nmax[k];
    dnspin_count += (nmax+1);
    if (nmax != -1) kshells_dn_.push_back({k,0,nmax});
  }
  if ((upspin_count+dnspin_count) != num_particle) {
    throw std::logic_error("* Fermisea::construct_groundstate:: consistency check-3 failed!");
  }

  // check
  /*
  for (int i=0; i<kshells_up_.size(); ++i) {
    std::cout<<std::setw(4)<<kshells_up_[i].k<<std::setw(4)<<kshells_up_[i].nmin
      <<std::setw(4)<<kshells_up_[i].nmax<<"\n";
  }
  getchar();
  for (int i=0; i<kshells_dn_.size(); ++i) {
    std::cout<<std::setw(4)<<kshells_dn_[i].k<<std::setw(4)<<kshells_dn_[i].nmin
      <<std::setw(4)<<kshells_dn_[i].nmax<<"\n";
  }
  getchar();
  */

  std::cout.copyfmt(state);
}

void Fermisea::get_sc_correlation(void)
{
  corr_aa_.setZero();
  corr_ab_.setZero();
  corr_fs_.setZero();
  for (int r=1; r<rmax_; ++r) {
    Vector3d R = R_list_[r];
    for (int i=0; i<kshells_up_.size(); ++i) {
      int k = kshells_up_[i].k;
      Vector3d kvec = blochbasis_.kvector(k);
      auto k_ab = std::exp(-ii()*kvec.dot(alpha_-beta_));
      // 1-body correlations
      for (int n=0; n<kblock_dim_; ++n) {
        corr_fs_(n,r) += std::exp(-ii()*kvec.dot(R));
      }
      for (int j=0; j<kshells_up_.size(); ++j) {
        int kp = kshells_up_[j].k;
        Vector3d kpvec = blochbasis_.kvector(kp);
        Vector3d kkplus = kvec+kpvec;
        Vector3d kkminus = kvec-kpvec;
        auto kkp = std::exp(ii()*kkplus.dot(R));
        auto kp_ab = std::exp(ii()*kpvec.dot(beta_-alpha_));
        auto term1 = 2*std::cos(kkminus.dot(alpha_));
        auto term2 = std::exp(ii()*(-kvec.dot(alpha_)+kpvec.dot(beta_)))
                    +std::exp(ii()*(kvec.dot(beta_)-kpvec.dot(alpha_)));
        for (int n=0; n<kblock_dim_; ++n) {
          if (n<=kshells_up_[i].nmax && n<=kshells_up_[j].nmax) {
            corr_aa_(n,r) += kkp*(2.0 + term1);
            corr_ab_(n,r) += kkp*(k_ab + term2 + kp_ab);
          }
        }
      }
    }
    for (int n=0; n<kblock_dim_; ++n) {
      corr_aa_(n,r) /= (2*num_kpoints_*num_kpoints_);
      corr_ab_(n,r) /= (2*num_kpoints_*num_kpoints_);
      corr_fs_(n,r) /= num_kpoints_;
    }
  }

  // print 
  std::ofstream fs("fs_sccorr.txt");
  fs << std::right<<std::scientific<<std::uppercase<<std::setprecision(6);
  for (int r=1; r<rmax_; ++r) {
    fs<<std::setw(6)<<r<< " ";
    for (int n=0; n<kblock_dim_; ++n) {
      //fs<<std::setw(15)<<std::real(corr_fs_(n,r)); 
      //fs<<std::setw(15)<<std::imag(corr_fs_(n,r)); 
      fs<<std::setw(15)<<std::real(corr_aa_(n,r)); 
      fs<<std::setw(15)<<std::imag(corr_aa_(n,r)); 
      fs<<std::setw(15)<<std::real(corr_ab_(n,r)); 
      fs<<std::setw(15)<<std::imag(corr_ab_(n,r)); 
    }
    fs<<"\n";
  }
  fs.close();
  std::cout << "sc_correlation done\n";
  //getchar();
}


} // end namespace var

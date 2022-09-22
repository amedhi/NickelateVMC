/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 23:06:41
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 11:32:25
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "strmatrix.h"
#include <boost/algorithm/string.hpp>
#include "./bcs_state.h"

namespace var {

BCS_State::BCS_State(const MF_Order::order_t& order, const MF_Order::pairing_t& pair_symm,
    const input::Parameters& inputs, const lattice::Lattice& lattice)
  : GroundState(order,pair_symm)
{
  init(inputs, lattice);
}

int BCS_State::init(const input::Parameters& inputs, const lattice::Lattice& lattice)
{
  name_ = "BCS";
  lattice_id_ = lattice.id();
  // sites & bonds
  num_sites_ = lattice.num_sites();
  num_bonds_ = lattice.num_bonds();
  // particle number
  set_particle_num(inputs);
  // infinity limit
  large_number_ = 5.0E+4;

  // Bloch basis
  blochbasis_.construct(lattice);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(lattice);

  // build MF Hamiltonian
  varparms_.clear();
  mf_model_.init(lattice);
  std::string name, path;
  double defval, lb, ub, dh;
  using namespace model;
  model::CouplingConstant cc;

  using order_t = MF_Order::order_t;
  using pairing_t = MF_Order::pairing_t;
  bool mf_model_finalized = false;
  int info;
  // read parameter
  strMatrix expr_mat;

  wf_analytical_form_ = inputs.set_value("wf_analytical_form",false,info);
  interband_pairing_ = inputs.set_value("interband_pairing",false,info);
//---------------------------------------------------------------------------
  if (lattice.id()==lattice::lattice_id::SQUARE) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    if (order()==order_t::SC && pair_symm()==pairing_t::SWAVE) {
      order_name_ = "SC-SWAVE";
      mf_model_.add_siteterm(name="pairing", cc="delta_sc", op::pair_create());
    }
    else if (order()==order_t::SC && pair_symm()==pairing_t::DWAVE) {
      order_name_ = "SC-DWAVE";
      cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    }
    else if (order()==order_t::SC && pair_symm()==pairing_t::EXTENDED_S) {
      order_name_ = "SC-Extended_S";
      cc = CouplingConstant({0, "delta_sc"}, {1, "delta_sc"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for the lattice");
    }
    // variational parameters
    defval = mf_model_.get_parameter_value("delta_sc");
    varparms_.add("delta_sc",defval,lb=1.0E-4,ub=6.0,dh=0.02);
    add_chemical_potential(inputs);
  }
  //---------------------------------------------------------------------------

  else if (lattice.id()==lattice::lattice_id::SQUARE_2SITE) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=0.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
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

    if (order()==order_t::SC && pair_symm()==pairing_t::DWAVE) {
      order_name_ = "SC-DWAVE";
      cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"}, {2,"0"}, {3,"0"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    }
    else if (order()==order_t::SC && pair_symm()==pairing_t::EXTENDED_S) {
      order_name_ = "SC-Extended_S";
      cc = CouplingConstant({0, "delta_sc"}, {1, "delta_sc"}, {2,"0"}, {3,"0"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for the lattice");
    }

    // variational parameters
    defval = mf_model_.get_parameter_value("delta_sc");
    varparms_.add("delta_sc", defval,lb=1.0E-3,ub=2.0,dh=0.02);
    defval = mf_model_.get_parameter_value("delta_af");
    varparms_.add("delta_af",defval,lb=0.0,ub=2.0,dh=0.02);
    add_chemical_potential(inputs);
  }

  else if (lattice.id()==lattice::lattice_id::SQUARE_STRIPE) {
    mf_model_.add_parameter(name="tx", defval=1.0, inputs);
    mf_model_.add_parameter(name="ty1", defval=1.0, inputs);
    mf_model_.add_parameter(name="ty2", defval=1.0, inputs);
    mf_model_.add_parameter(name="lambda", defval=8, inputs);
    // variational parameter
    mf_model_.add_parameter(name="dSC_x", defval=1.0, inputs);
    mf_model_.add_parameter(name="dSC_y1", defval=1.0, inputs);
    mf_model_.add_parameter(name="dSC_y2", defval=1.0, inputs);
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

    if (order()==order_t::SC && pair_symm()==pairing_t::SC_CDW_SDW) {
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

      // pairing term
      cc.create(48);
      // x-bonds
      for (int n=0; n<16; ++n) {
        double ampl;
        if (lambda > 0) ampl = std::abs(std::cos(0.5*Q*(n+0.5)));
        else ampl = 1.0;
        strval << ampl << "*dSC_x";
        cc.add_type(n, strval.str());
        //std::cout << "ampl["<<n<<"] = "<< strval.str() << "\n";
        strval.str(std::string()); // clear
      }
      // y-bonds (intracell)
      for (int n=16; n<32; ++n) {
        int m = n-16;
        double ampl;
        if (lambda > 0) ampl = -std::abs(std::cos(0.5*Q*(m+0.5)));
        else ampl = -1.0;
        strval << ampl << "*dSC_y1";
        cc.add_type(n, strval.str());
        //std::cout << "ampl["<<n<<"] = "<< strval.str() << "\n";
        strval.str(std::string()); // clear
      }
      // y-bonds (intercell)
      for (int n=32; n<48; ++n) {
        int m = n-32;
        double ampl;
        if (lambda>0) ampl = -std::abs(std::cos(0.5*Q*(m+0.5)));
        else ampl = -1.0;
        strval << ampl << "*dSC_y2";
        cc.add_type(n, strval.str());
        //std::cout << "ampl["<<n<<"] = "<< strval.str() << "\n";
        strval.str(std::string()); // clear
      }
      // pairing term
      mf_model_.add_bondterm(name="bond_singlet", cc, op::pair_create());

      // variational parameters
      defval = mf_model_.get_parameter_value("dSC_x");
      varparms_.add("dSC_x", defval,lb=1.0E-3,ub=2.0,dh=0.01);
      defval = mf_model_.get_parameter_value("dSC_y1");
      varparms_.add("dSC_y1", defval,lb=1.0E-3,ub=2.0,dh=0.01);
      defval = mf_model_.get_parameter_value("dSC_y2");
      varparms_.add("dSC_y2", defval,lb=1.0E-3,ub=2.0,dh=0.01);
      defval = mf_model_.get_parameter_value("dCDW");
      varparms_.add("dCDW", defval,lb=0,ub=2.0,dh=0.01);
      defval = mf_model_.get_parameter_value("dSDW");
      varparms_.add("dSDW", defval,lb=0,ub=2.0,dh=0.01);

      add_chemical_potential(inputs);
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for the lattice");
    }
  }


//---------------------------------------------------------------------------
  else if (lattice.id()==lattice::lattice_id::HONEYCOMB) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
    mf_model_.add_parameter(name="theta1", defval=0.0, inputs);
    mf_model_.add_parameter(name="theta2", defval=0.0, inputs);
    // bond operators
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
    // site operators
    //mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
    // pairing term
    cc = CouplingConstant({0,"delta_sc"}, {1,"delta_sc*exp(i*theta1)"}, 
        {2,"delta_sc*exp(i*theta2)"});
    mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
    // variational parameters
    defval = mf_model_.get_parameter_value("delta_sc");
    varparms_.add("delta_sc", defval, lb=1.0E-3, ub=2.0, dh=0.01);
    defval = mf_model_.get_parameter_value("theta1");
    varparms_.add("theta1", defval, lb=0.0, ub=two_pi(), dh=0.2);
    defval = mf_model_.get_parameter_value("theta2");
    varparms_.add("theta2", defval, lb=0.0, ub=two_pi(), dh=0.2);

    add_chemical_potential(inputs);
  }

//---------------------------------------------------------------------------
  else if (lattice.id()==lattice::lattice_id::NICKELATE_2L) {
    bool mu_default = inputs.set_value("mu_default", false);
    mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=0.0, inputs);
    mf_model_.add_parameter(name="th", defval=0.0, inputs);
    mf_model_.add_parameter(name="delta_N", defval=1.0, inputs);
    mf_model_.add_parameter(name="delta_R", defval=1.0, inputs);
    mf_model_.add_parameter(name="mu_N", defval=0.0, inputs, info);
    if (mu_default && info==0) std::cout << " >> alert: conflicting option for `mu`\n";
    mf_model_.add_parameter(name="mu_R", defval=0.0, inputs, info);
    if (mu_default && info==0) std::cout << " >> alert: conflicting option for `mu`\n";

    // site operators
    cc.create(2);
    cc.add_type(0, "-mu_N");
    cc.add_type(1, "e_R-mu_R");
    //cc.add_type(1, "e_R-mu_N");
    mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());

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

    // pairing term
    if (order()==order_t::SC && pair_symm()==pairing_t::SWAVE) {
      order_name_ = "SC-SWAVE";
      cc.create(2);
      cc.add_type(0, "delta_N");
      cc.add_type(1, "delta_R");
      //cc.add_type(1, "0");
      mf_model_.add_siteterm(name="singlet", cc, op::pair_create());
    }

    else if (order()==order_t::SC && pair_symm()==pairing_t::DWAVE) {
      order_name_ = "SC-DWAVE";
      cc.create(7);
      cc.add_type(0, "delta_N");
      cc.add_type(1, "-delta_N");
      cc.add_type(2, "0");
      cc.add_type(3, "delta_R");
      cc.add_type(4, "-delta_R");
      cc.add_type(5, "0");
      cc.add_type(6, "0");
      mf_model_.add_bondterm(name="bond_singlet", cc, op::pair_create());
    }

    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for the lattice");
    }
    // variational parameters
    defval = mf_model_.get_parameter_value("delta_N");
    varparms_.add("delta_N", defval,lb=1.0E-3,ub=6.0,dh=0.02);
    defval = mf_model_.get_parameter_value("delta_R");
    varparms_.add("delta_R",defval,lb=0.0,ub=6.0,dh=0.02);

    // chemical potential
    noninteracting_mu_ = false;
    if (inputs.set_value("mu_default", false)) {
      // starting with non-interacting mu
      mf_model_.finalize(lattice);
      mf_model_.update_site_parameter("mu_N", 0.0);
      mf_model_.update_site_parameter("mu_R", 0.0);
      double mu_0 = get_noninteracting_mu();
      //std::cout << "mu = " << mu_0 << "\n"; getchar();
      mf_model_.update_site_parameter("mu_N", mu_0);
      mf_model_.update_site_parameter("mu_R", mu_0);
      mf_model_finalized = true;
    }
    defval = mf_model_.get_parameter_value("mu_N");
    varparms_.add("mu_N",defval,lb=defval-10.0,ub=defval+10.0,dh=0.1);
    defval = mf_model_.get_parameter_value("mu_R");
    varparms_.add("mu_R",defval,lb=defval-10.0,ub=defval+10.0,dh=0.1);
  }
//---------------------------------------------------------------------------
  else {
    throw std::range_error("BCS_State::BCS_State: undefined for the lattice");
  }

  // finalize MF Hamiltonian
  num_varparms_ = varparms_.size();
  if (!mf_model_finalized) {
    mf_model_.finalize(lattice);
  }

  //std::cout<< "mu_0 = " <<  get_noninteracting_mu() << "\n";

  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  delta_k_.resize(kblock_dim_,kblock_dim_);
  dphi_k_.resize(kblock_dim_,kblock_dim_);
  phi_k_.resize(num_kpoints_);
  work_k_.resize(num_kpoints_);
  bdg_mat_.resize(2*kblock_dim_,2*kblock_dim_);
  uk_.resize(kblock_dim_,kblock_dim_);
  vk_.resize(kblock_dim_,kblock_dim_);
  ubk_.resize(kblock_dim_,kblock_dim_);
  vbk_.resize(kblock_dim_,kblock_dim_);
  for (int k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
    work_k_[k].resize(kblock_dim_,kblock_dim_);
  } 
  return 0;
}

void BCS_State::update(const lattice::Lattice& lattice)
{
  // update for change in lattice BC (same structure & size)
  // bloch basis
  blochbasis_.construct(lattice);
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(lattice);
}

void BCS_State::add_chemical_potential(const input::Parameters& inputs)
{
  // default handling
  int info;
  double defval, lb, ub, dh;
  std::string name;
  model::CouplingConstant cc;
  using namespace model;
  mf_model_.add_parameter(name="mu", defval=0.0, inputs, info);
  mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
  if (info == 0) noninteracting_mu_ = false;
  else noninteracting_mu_ = true;
  if (inputs.set_value("mu_variational", false, info)) {
    varparms_.add("mu",defval=0.0,lb=-5.0,ub=+5.0,dh=0.1);
  }
}

std::string BCS_State::info_str(void) const
{
  std::ostringstream info;
  if (kblock_dim_>1) {
    if (interband_pairing_) {
      info << "# Ground State: '"<<name_<<" ("<<order_name_<<") (INTERBAND)'\n";
    }
    else {
      info << "# Ground State: '"<<name_<<" ("<<order_name_<<") (INTRABAND)'\n";
    }
  }
  else {
    info << "# Ground State: '"<<name_<<" ("<<order_name_<<")'\n";
  }
  info << "# Hole doping = "<<hole_doping()<<"\n";
  info << "# Particles = "<< num_upspins()+num_dnspins();
  info << " (Nup = "<<num_upspins()<<", Ndn="<<num_dnspins()<<")\n";
  info.precision(6);
  info.setf(std::ios_base::fixed);
  if (noninteracting_mu_)
    info << "# mu = non-interacting value\n";
  //else info << "# mu = "<<mf_model_.get_parameter_value("mu")<<"\n";
  return info.str();
}

void BCS_State::update(const input::Parameters& inputs)
{
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
  // update MF model
  mf_model_.update(inputs);
  // update variational parameters
  for (auto& p : varparms_) 
    p.change_value(mf_model_.get_parameter_value(p.name()));
  // chemical potential
  if (noninteracting_mu_) {
    /*int nowarn;
    if (inputs.set_value("custom_mu", false, nowarn)) {
      mf_model_.update_site_parameter("mu_N", 0.0);
      mf_model_.update_site_parameter("mu_R", 0.0);
      double mu_0 = get_noninteracting_mu();
      std::cout << "mu_0 = " << mu_0 << "\n";
      mf_model_.update_site_parameter("mu_N", mu_0);
      mf_model_.update_site_parameter("mu_R", mu_0);
    }
    else {
      // the next line is 'really' needed 
      mf_model_.update_site_parameter("mu", 0.0);
      double mu_0 = get_noninteracting_mu();
      //std::cout << "mu = " << mu_0 << "\n";
      mf_model_.update_site_parameter("mu", mu_0);
    }
    */
    // the next line is 'really' needed 
    mf_model_.update_site_parameter("mu", 0.0);
    double mu_0 = get_noninteracting_mu();
    std::cout << "mu = " << mu_0 << "\n";
    mf_model_.update_site_parameter("mu", mu_0);
  }
  // check MF energy
  //std::cout << " MF energy = " << get_mf_energy() << "\n"; 
}

void BCS_State::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
  // update from variational parameters
  int i = 0;
  for (auto& p : varparms_) {
    auto x = pvector[start_pos+i];
    p.change_value(x);
    mf_model_.update_parameter(p.name(), x);
    //std::cout << p.name() << " = " << x << "\n";
    ++i;
  }
  mf_model_.update_terms();
}

void BCS_State::get_wf_amplitudes(Matrix& psi) 
{
  // k-space pair amplitudes
  if (wf_analytical_form_) {
    if (lattice_id_==lattice::lattice_id::NICKELATE_2L) {
      analytical_amplitudes_NICKELATE2L(phi_k_);
    }
    else {
      throw std::range_error("BCS_State: analytical wf not defined for this lattice");
    }
  }
  else {
    if (kblock_dim_==1) {
      get_pair_amplitudes_oneband(phi_k_);
    }
    else {
      if (interband_pairing_) {
        get_pair_amplitudes_interband(phi_k_);
      }
      else {
        get_pair_amplitudes_intraband(phi_k_);
      }
    }
  }

  // 'lattice space' pair amplitudes
  get_pair_amplitudes_sitebasis(phi_k_, psi);
}

void BCS_State::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  // k-space pair amplitudes
  if (wf_analytical_form_) {
    if (lattice_id_==lattice::lattice_id::NICKELATE_2L) {
      analytical_gradient_NICKELATE2L(psi_gradient);
    }
    else {
      throw std::range_error("BCS_State: analytical wf not defined for this lattice");
    }
    return;
  }
  //std::cout << "Numerical gradient\n";
  // numerical gradient
  int i=0; 
  for (const auto& p : varparms_) {
    double h = p.diff_h();
    double x = p.value();
    double x_fwd = x+h; 
    if (x_fwd > p.ubound()) x_fwd = p.ubound();
    mf_model_.update_parameter(p.name(), x_fwd);
    mf_model_.update_terms();
    if (kblock_dim_==1) get_pair_amplitudes_oneband(phi_k_);
    else {
      if (interband_pairing_) {
        get_pair_amplitudes_interband(phi_k_);
      }
      else {
        get_pair_amplitudes_intraband(phi_k_);
      }
    } 
    double x_bwd = x-h; 
    if (x_bwd < p.lbound()) x_bwd = p.lbound();
    mf_model_.update_parameter(p.name(), x_bwd);
    mf_model_.update_terms();
    if (kblock_dim_==1) get_pair_amplitudes_oneband(work_k_);
    else {
      if (interband_pairing_) {
        get_pair_amplitudes_interband(work_k_);
      }
      else {
        get_pair_amplitudes_intraband(work_k_);
      }
    } 
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
    ++i;
  }
}

void BCS_State::get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
  // psi = FTU_ * PHI_K * conjugate(transpose(FTU_))
  // PHI_K is block diagonal (k-th block is phi_k) 
  int p = 0;
  for (int i=0; i<num_kpoints_; ++i) {
    int q = 0;
    for (int j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (int k=0; k<num_kpoints_; ++k) {
        work_ += FTU_(i,k) * phi_k[k] * std::conj(FTU_(j,k));
      }
      // copy transformed block
      //psi.block(p,q,kblock_dim_,kblock_dim_) = 
      for (int m=0; m<kblock_dim_; ++m) {
        for (int n=0; n<kblock_dim_; ++n) 
          psi(p+m,q+n) = ampl_part(work_(m,n));
      }
      q += kblock_dim_;
    }
    p += kblock_dim_;
  }
  /*
  for (int i=0; i<num_sites_; ++i) {
    for (int j=0; j<num_sites_; ++j) {
      std::cout << "psi["<<i<<","<<j<<"] = "<<psi(i,j)<<"\n";
      getchar();
    }
  }
  */
}

void BCS_State::get_pair_amplitudes_oneband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for one band system 
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    double ek = std::real(mf_model_.quadratic_spinup_block()(0,0)); 
    auto delta_k = mf_model_.pairing_part()(0,0);
    //std::cout<<"delta_k="<<delta_k<<"\n"; getchar();
    mf_model_.construct_kspace_block(-kvec);
    ek += std::real(mf_model_.quadratic_spinup_block()(0,0));;
    delta_k += mf_model_.pairing_part()(0,0);
    delta_k *= 0.5;
    //----------------------------------
    //std::cout << "** Hack at BCS_State\n";
    //ek = -4.0*(std::cos(kvec[0])+std::cos(kvec[1]));
    //delta_k = 0.5 * (std::cos(kvec[0])-std::cos(kvec[1])); 
    //----------------------------------
    double deltak_sq = std::norm(delta_k);
    double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
    if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
      phi_k[k](0,0) = large_number_ * std::exp(ii()*std::arg(delta_k));
    }
    else {
      phi_k[k](0,0) = 2.0*delta_k/ek_plus_Ek;
    }
    /*
    std::cout << "k = " << k << "\n";
    std::cout << "kvec = " << kvec.transpose() << "\n";
    std::cout << "ek = " << 0.5*ek << "\n";
    std::cout << "delta_k = " << delta_k << "\n";
    std::cout << "phi_k = " << phi_k[k](0,0) << "\n";
    getchar();
    */
  }
}

void BCS_State::get_pair_amplitudes_intraband(std::vector<ComplexMatrix>& phi_k)
{
  //std::cout << "get_pair_amplitudes_intraband\n";
  // BCS pair amplitudes for multi-band system (INTRABAND pairing only)
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    //std::cout << es_k_up.eigenvalues().transpose() << "\n";
    //getchar();

    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //std::cout << "kvec = " << kvec.transpose() << "\n";
    //std::cout << "delta_k =\n" << delta_k_ << "\n\n";

    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_dn.compute(mf_model_.quadratic_spindn_block());

    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    //std::cout << "delta_ck =\n" << work_ << "\n\n";

    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_dn.eigenvectors().conjugate();
    //std::cout << "Uk_up =\n" << es_k_up.eigenvectors() << "\n\n";
    //std::cout << "delta_k_ =\n" << delta_k_.real() << "\n"; getchar();

    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    dphi_k_.setZero();
    for (int i=0; i<kblock_dim_; ++i) {
      double ek = es_k_up.eigenvalues()[i] + es_minusk_dn.eigenvalues()[i];
      double deltak_sq = std::norm(delta_k_(i,i));
      double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
      /*
      if (ek<0.0) dphi_k_(i,i) = 1.0;
      else dphi_k_(i,i) = 0.0;
      */
      if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
        dphi_k_(i,i) = large_number_ * std::exp(ii()*std::arg(delta_k_(i,i)));
        //std::cout << "n = " << i << "\n";
        //std::cout << "kvec = " << kvec.transpose() << "\n";
        //std::cout << "deltak_sq = "<<deltak_sq<<"\n";
        //std::cout << "ek = "<<ek<<"\n";
        //std::cout << "ek_plus_Ek = "<<ek_plus_Ek<<"\n";
        std::cout << ">> alert! BCS_State: singularilty in 'phi_k', replaced by 'large number'\n";
        //getchar();
      }
      else {
        dphi_k_(i,i) = 2.0*delta_k_(i,i)/ek_plus_Ek;
      }
      /*
      std::cout << "ek = " << ek << "\n";
      std::cout << "delta_sq = " << deltak_sq << "\n";
      std::cout << "phi_k = " << dphi_k_(i,i) << "\n";
      getchar();
      */
    }
    // bcs ampitudes in original basis 
    for (int i=0; i<kblock_dim_; ++i) 
      work_.col(i) = es_k_up.eigenvectors().col(i) * dphi_k_(i,i);
    phi_k[k] = work_ * es_minusk_dn.eigenvectors().transpose();
  } 
}

void BCS_State::get_pair_amplitudes_interband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for multi-band system (Including INTERBAND pairing)
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    //std::cout << es_k_up.eigenvalues().transpose() << "\n";
    //getchar();

    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //std::cout << "kvec = " << kvec.transpose() << "\n";
    //std::cout << "delta_k =\n" << delta_k_ << "\n\n";

    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_dn.compute(mf_model_.quadratic_spindn_block());

    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    //std::cout << "delta_ck =\n" << work_ << "\n\n";

    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_dn.eigenvectors().conjugate();
    //std::cout << "Uk_up =\n" << es_k_up.eigenvectors() << "\n\n";
    //std::cout << "delta_k_ =\n" << delta_k_.real() << "\n"; getchar();

    // BCS ampitudes in rotated basis 
    int dim = kblock_dim_;
    bdg_mat_.setZero();
    // Construct the BDG matrix
    for (int i=0; i<dim; ++i) {
      bdg_mat_(i,i) = es_k_up.eigenvalues()[i];
      bdg_mat_(dim+i,dim+i) = -es_minusk_dn.eigenvalues()[i];
    }
    bdg_mat_.block(0,dim, dim,dim) = -delta_k_;
    bdg_mat_.block(dim,0, dim,dim) = -delta_k_.adjoint();
    //std::cout << "BDG matrix\n";
    //std::cout << bdg_mat_.real() << "\n\n";
    // Diagonalize BDG matrix
    es_bdg_.compute(bdg_mat_);
    ComplexMatrix Udag = es_bdg_.eigenvectors().adjoint();
    //std::cout << "Udag =\n" << Udag.real() << "\n\n";

    // uk & vk matrices (since +ve energy Bogoliubov operators in 2nd half)
    uk_ = Udag.block(dim,0,dim,dim);
    vk_ = Udag.block(dim,dim,dim,dim);

    // phi_k = uk^{-1} * vk
    // SVD decompose: uk = USV^\dag => (uk)^{-1} = VS^{-1}U^{\dag} 
    Eigen::JacobiSVD<Eigen::MatrixXcd> svd(uk_, Eigen::ComputeThinU | Eigen::ComputeThinV); 
    work_ = svd.matrixU().adjoint() * vk_;
    bool singular = false;
    for (int i=0; i<dim; ++i) {
      double dval = svd.singularValues()(i);
      if (std::abs(dval) < 1.0E-12) {
        //std::cout << uk_.real() << "\n";
        //std::cout << vk_.row(i).real() << "\n";
        //std::cout << "s["<<i<<"] = " << dval << "\n";

        /*
        std::cout << ">> alert! BCS_State: uk matrix singular for k = (";
        std::cout << kvec[0]<<", "<<kvec[1]<<", "<<kvec[2]<<")\n";
        */

        //getchar();
        singular = true;
        break;
      }
      work_.row(i) /= dval;
    } 
    if (!singular) {
      dphi_k_ = svd.matrixV() * work_;
    }
    else {
      // INTRABAND pairing only for this k
      dphi_k_.setZero();
      for (int i=0; i<kblock_dim_; ++i) {
        double ek = es_k_up.eigenvalues()[i] + es_minusk_dn.eigenvalues()[i];
        double deltak_sq = std::norm(delta_k_(i,i));
        double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
        if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
          dphi_k_(i,i) = large_number_ * std::exp(ii()*std::arg(delta_k_(i,i)));
          //std::cout << ">> alert! BCS_State: singularilty in 'phi_k', replaced by 'large number'\n";
          /*
          std::cout << "deltak_sq = "<<deltak_sq<<"\n";
          std::cout << "ek = "<<ek<<"\n";
          std::cout << "ek_plus_Ek = "<<ek_plus_Ek<<"\n";
          getchar();
          */
        }
        else {
          dphi_k_(i,i) = 2.0*delta_k_(i,i)/ek_plus_Ek;
        }
      }
    }
    // phi_k for 'ck'-operator
    phi_k[k] = es_k_up.eigenvectors()*dphi_k_*es_minusk_dn.eigenvectors().transpose();
  } 
}

double BCS_State::get_mf_energy(void)
{
  double mf_energy = 0.0;
  double delta = mf_model_.get_parameter_value("delta_sc");
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    /*
    double ek = -2.0*(std::cos(kvec[0])+std::cos(kvec[1]));
    double deltak = delta * (std::cos(kvec[0])-std::cos(kvec[1]));
    double Ek = std::sqrt(ek*ek + deltak*deltak);
    mf_energy += (ek - Ek);
    */
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_dn.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    //std::cout << work_ << "\n"; getchar();
    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_dn.eigenvectors().conjugate();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    for (unsigned i=0; i<kblock_dim_; ++i) {
      double ek = es_k_up.eigenvalues()[i];
      double deltak_sq = std::norm(delta_k_(i,i));
      double Ek = std::sqrt(ek*ek + deltak_sq);
      //std::cout << ek << " " << Ek << "\n"; getchar();
      mf_energy += (ek - Ek);
    }
  }
  // constant term
  double J = 0.35;  // tJ model
  mf_energy += num_bonds_*delta*delta/(2.0*J);  
  //std::cout << delta << " ";
  return mf_energy/num_sites_;
}

void BCS_State::analytical_amplitudes_NICKELATE2L(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes, lattice NICKELATE_2L (INTRABAND pairing only)
  if (pair_symm()!=pairing_t::DWAVE) {
    throw std::range_error("BCS_State: analytical wf not defined for this pairing symmetry");
  }
  // read parameters
  double t = mf_model_.get_parameter_value("t");
  double tp = mf_model_.get_parameter_value("tp");
  double th = mf_model_.get_parameter_value("th");
  double eR = mf_model_.get_parameter_value("e_R");
  double DeltaN = mf_model_.get_parameter_value("delta_N");
  double DeltaR = mf_model_.get_parameter_value("delta_R");
  double muN = mf_model_.get_parameter_value("mu_N");
  double muR = mf_model_.get_parameter_value("mu_R");

  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    double kx = kvec[0];
    double ky = kvec[1];
    double ek = -2.0*t*(std::cos(kx)+std::cos(ky))-4.0*tp*std::cos(kx)*std::cos(ky); 
    double eps1k = ek - muN;
    double eps2k = ek + eR - muR;
    double e1k = eps2k + eps1k;
    double e2k = eps2k - eps1k;
    double exk = std::sqrt(e2k*e2k + 4.0*th*th);
    double ak = (e2k + exk)/(2.0*th);
    double xk = std::sqrt(1.0+ak*ak);
    double w1k = ak/xk;
    double w2k = 1.0/xk;
    double DeltaN_k = DeltaN*(std::cos(kx)-std::cos(ky));
    double DeltaR_k = DeltaR*(std::cos(kx)-std::cos(ky));
    double Delta1k = w1k*w1k*DeltaN_k + w2k*w2k*DeltaR_k;
    double Delta2k = w1k*w1k*DeltaR_k + w2k*w2k*DeltaN_k;
    double E1k = 0.5*(e1k - exk);
    double E2k = 0.5*(e1k + exk);
    // phi_k(1) & phi_k(2)
    double Delta1k_sq = Delta1k*Delta1k;
    double Delta2k_sq = Delta2k*Delta2k;
    double Eqp1k = std::sqrt(E1k*E1k + Delta1k_sq);
    double Eqp2k = std::sqrt(E2k*E2k + Delta2k_sq);
    dphi_k_.setZero();
    if (std::sqrt(Delta1k_sq)<1.0E-12 && E1k<0.0) {
      dphi_k_(0,0) = large_number_;
    }
    else {
      dphi_k_(0,0) = Delta1k/(E1k + Eqp1k);
    }
    if (std::sqrt(Delta2k_sq)<1.0E-12 && E2k<0.0) {
      dphi_k_(1,1) = large_number_;
    }
    else {
      dphi_k_(1,1) = Delta2k/(E2k + Eqp2k);
    }

    // bcs ampitudes in original basis 
    ComplexMatrix Uk(2,2);
    Uk(0,0) = w1k; Uk(0,1) = -w2k;
    Uk(1,0) = w2k; Uk(1,1) =  w1k;
    phi_k[k] = Uk * dphi_k_ * Uk.transpose();
  }
}

void  BCS_State::analytical_gradient_NICKELATE2L(std::vector<Matrix>& psi_gradient)
{
  // BCS pair amplitudes, lattice NICKELATE_2L (INTRABAND pairing only)
  if (pair_symm()!=pairing_t::DWAVE) {
    throw std::range_error("BCS_State: analytical gradient not defined for this pairing symmetry");
  }
  std::cout << "Analytical gradient\n";
  // read parameters
  double t = mf_model_.get_parameter_value("t");
  double tp = mf_model_.get_parameter_value("tp");
  double th = mf_model_.get_parameter_value("th");
  double eR = mf_model_.get_parameter_value("e_R");
  double DeltaN = mf_model_.get_parameter_value("delta_N");
  double DeltaR = mf_model_.get_parameter_value("delta_R");
  double muN = mf_model_.get_parameter_value("mu_N");
  double muR = mf_model_.get_parameter_value("mu_R");

  // storage
  std::vector<ComplexMatrix> Dphik_deltaN(num_kpoints_);
  std::vector<ComplexMatrix> Dphik_deltaR(num_kpoints_);
  std::vector<ComplexMatrix> Dphik_muN(num_kpoints_);
  std::vector<ComplexMatrix> Dphik_muR(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    Dphik_deltaN[k].resize(kblock_dim_,kblock_dim_);
    Dphik_deltaR[k].resize(kblock_dim_,kblock_dim_); 
    Dphik_muN[k].resize(kblock_dim_,kblock_dim_);
    Dphik_muR[k].resize(kblock_dim_,kblock_dim_);
  } 

  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    double kx = kvec[0];
    double ky = kvec[1];
    double ek = -2.0*t*(std::cos(kx)+std::cos(ky))-4.0*tp*std::cos(kx)*std::cos(ky); 
    double eps1k = ek - muN;
    double eps2k = ek + eR - muR;
    double e1k = eps2k + eps1k;
    double e2k = eps2k - eps1k;
    double exk = std::sqrt(e2k*e2k + 4.0*th*th);
    double ak = (e2k + exk)/(2.0*th);
    double xk = std::sqrt(1.0+ak*ak);
    double w1k = ak/xk;
    double w2k = 1.0/xk;
    double coskx_minus_cosky = (std::cos(kx)-std::cos(ky));
    double DeltaN_k = DeltaN * coskx_minus_cosky;
    double DeltaR_k = DeltaR * coskx_minus_cosky;
    double Delta1k = w1k*w1k*DeltaN_k + w2k*w2k*DeltaR_k;
    double Delta2k = w1k*w1k*DeltaR_k + w2k*w2k*DeltaN_k;
    double E1k = 0.5*(e1k - exk);
    double E2k = 0.5*(e1k + exk);
    // phi_k(1) & phi_k(2)
    double Delta1k_sq = Delta1k*Delta1k;
    double Delta2k_sq = Delta2k*Delta2k;
    double Eqp1k = std::sqrt(E1k*E1k + Delta1k_sq);
    double Eqp2k = std::sqrt(E2k*E2k + Delta2k_sq);
    dphi_k_.setZero();
    if (std::sqrt(Delta1k_sq)<1.0E-12 && E1k<0.0) {
      dphi_k_(0,0) = large_number_;
    }
    else {
      dphi_k_(0,0) = Delta1k/(E1k + Eqp1k);
    }
    if (std::sqrt(Delta2k_sq)<1.0E-12 && E2k<0.0) {
      dphi_k_(1,1) = large_number_;
    }
    else {
      dphi_k_(1,1) = Delta2k/(E2k + Eqp2k);
    }

    // simplifying
    double phi1k = std::real(dphi_k_(0,0));
    double phi2k = std::real(dphi_k_(1,1));

    double Dphi_delta1 = phi1k*E1k/(Delta1k*Eqp1k);
    double Dphi_delta2 = phi2k*E2k/(Delta2k*Eqp2k);
    double w1k_sq = w1k*w1k;
    double w2k_sq = w2k*w2k;
    double Ddelta1_deltaN = w1k_sq * coskx_minus_cosky;
    double Ddelta2_deltaN = w2k_sq * coskx_minus_cosky;
    double Ddelta1_deltaR = Ddelta2_deltaN;
    double Ddelta2_deltaR = Ddelta1_deltaN;

    // For ampitudes in original basis 
    ComplexMatrix Uk(2,2);
    Uk(0,0) = w1k; Uk(0,1) = -w2k;
    Uk(1,0) = w2k; Uk(1,1) =  w1k;

    // Derivatives wrt 'Delta_N'
    dphi_k_(0,0) = Dphi_delta1 * Ddelta1_deltaN;
    dphi_k_(1,1) = Dphi_delta2 * Ddelta2_deltaN;
    Dphik_deltaN[k] = Uk * dphi_k_ * Uk.transpose();

    // Derivatives wrt 'Delta_R'
    dphi_k_(0,0) = Dphi_delta1 * Ddelta1_deltaR;
    dphi_k_(1,1) = Dphi_delta2 * Ddelta2_deltaR;
    Dphik_deltaR[k] = Uk * dphi_k_ * Uk.transpose();

    // Derivatives wrt 'mu_N'
    double Dphi_E1 = -phi1k/Eqp1k;
    double Dphi_E2 = -phi2k/Eqp2k;
    double DE1_muN = -0.5*(1.0+e2k/exk);
    double DE2_muN = -0.5*(1.0-e2k/exk);

    double Dw1_a = w2k/(1.0+ak*ak); 
    double Dw2_a = -ak*w2k/(1.0+ak*ak); 
    double Da_muN = ak/exk;
    double Dw1_muN = Dw1_a * Da_muN;
    double Dw2_muN = Dw2_a * Da_muN;
    double Ddelta1_muN = 2.0*w1k*DeltaN_k*Dw1_muN + 2.0*w2k*DeltaR_k*Dw2_muN;
    double Ddelta2_muN = 2.0*w2k*DeltaN_k*Dw2_muN + 2.0*w1k*DeltaR_k*Dw1_muN;

    dphi_k_(0,0) = Dphi_delta1*Ddelta1_muN + Dphi_E1*DE1_muN;
    dphi_k_(1,1) = Dphi_delta2*Ddelta2_muN + Dphi_E2*DE2_muN;
    Dphik_muN[k] = Uk * dphi_k_ * Uk.transpose();

    // Derivatives wrt 'mu_R'
    double DE1_muR = -0.5*(1.0-e2k/exk);
    double DE2_muR = -0.5*(1.0+e2k/exk);

    double Da_muR = 1.0/(ak*exk);
    double Dw1_muR = Dw1_a * Da_muR;
    double Dw2_muR = Dw2_a * Da_muR;
    double Ddelta1_muR = 2.0*w1k*DeltaN_k*Dw1_muR + 2.0*w2k*DeltaR_k*Dw2_muR;
    double Ddelta2_muR = 2.0*w2k*DeltaN_k*Dw2_muR + 2.0*w1k*DeltaR_k*Dw1_muR;

    dphi_k_(0,0) = Dphi_delta1*Ddelta1_muR + Dphi_E1*DE1_muR;
    dphi_k_(1,1) = Dphi_delta2*Ddelta2_muR + Dphi_E2*DE2_muR;
    Dphik_muR[k] = Uk * dphi_k_ * Uk.transpose();
  }

  // for all the variational parameters
  for (int p=0; p<varparms_.size(); ++p) {
    if (varparms_[p].name() == "delta_N") {
      get_pair_amplitudes_sitebasis(Dphik_deltaN, psi_gradient[p]);
    }
    else if (varparms_[p].name() == "delta_R") {
      get_pair_amplitudes_sitebasis(Dphik_deltaR, psi_gradient[p]);
    }
    else if (varparms_[p].name() == "mu_N") {
      get_pair_amplitudes_sitebasis(Dphik_muN, psi_gradient[p]);
    }
    else if (varparms_[p].name() == "mu_R") {
      get_pair_amplitudes_sitebasis(Dphik_muR, psi_gradient[p]);
    }
    else {
      throw std::range_error("BCS_State::analytical_gradient:: internal error");
    }
  }
}



} // end namespace var

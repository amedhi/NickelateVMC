/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:17:49
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-05 11:47:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <boost/algorithm/string.hpp>
#include "./projector.h"

namespace var {


/*----------------------------------------------------------------
* Gutzwiller Projection
*----------------------------------------------------------------*/
int GW_Projector::init(const lattice::Lattice& lattice, const input::Parameters& inputs,
  VariationalParms& vparms)
{
  // check if present
  int info;
  is_present_ = inputs.set_value("gw_projector",false,info);
  if (!is_present_) {
    default_case_ = false;
    return 0;
  }

  // initialize  
  num_sites_ = lattice.num_sites();
  num_site_types_ = lattice.num_site_types();
  basis_projection_.resize(num_site_types_);
  site_projection_.resize(num_sites_);
  for (auto& p : basis_projection_) p = pjn_t::NONE;
  for (auto& p : site_projection_) p = pjn_t::NONE;
  site_typeid_.resize(num_sites_);
  for (int i=0; i<num_sites_; ++i) {
    site_typeid_[i] = lattice.site(i).type();
  }
  gw_ratio_.resize(num_sites_,5);
  gw_ratio_.setConstant(1.0);


  // special case for 'NICKELATE_2B'
  if (lattice.id()==lattice::lattice_id::NICKELATE_2B) {
    // CAUTION: Assumption is that all the 'atom's are of the same type.
    nickelate_2B_ = true;
    default_case_ = false;
    uniform_projection_ = false;
    assert(num_site_types_==2);
    for (auto& p : basis_projection_) p = pjn_t::DOUBLON;
    for (auto& p : site_projection_) p = pjn_t::DOUBLON;

    // read gw factors
    gw_factor_.resize(3);
    gfactor_names_.resize(3);
    gfactor_names_[0] = "gfactor_ss";
    gfactor_names_[1] = "gfactor_dd";
    gfactor_names_[2] = "gfactor_sd";
    //gfactor_names_[3] = "gfactor_uu";
    for (int n=0; n<3; ++n) {
      gw_factor_[n] = inputs.set_value(gfactor_names_[n], 1.0);
      if (gw_factor_[n] < gw_cutoff()) {
        throw std::range_error("GW_Projector::init: 'gw_factor' value out of range");
      }
      // add to variational parameter list
      vparms.add(gfactor_names_[n], gw_factor_[n], gw_cutoff(), 10.0);
    }
    // store the 'factors' in matrix form
    num_orbitals_ = 2;
    gfactor_ud_.resize(2,2); 
    gfactor_uu_.resize(2,2); 
    gfactor_ud_(0,0) = gw_factor_[0];
    gfactor_ud_(1,1) = gw_factor_[1];
    gfactor_ud_(0,1) = gw_factor_[2];
    gfactor_ud_(1,0) = gw_factor_[2];
    //gfactor_uu_.setConstant(1.0);
    //gfactor_uu_(0,1) = gw_factor_[3];
    //gfactor_uu_(1,0) = gw_factor_[3];

    return gw_factor_.size();
  }


  // check if uniform or not
  if (num_site_types_ == 1) {
    uniform_projection_ = true;
  }
  else {
    // if 'gw_factor' not specified for every site types, then uniform
    uniform_projection_ = false;
    for (int n=0; n<num_site_types_; ++n) {
      std::string pname = "gw_factor"+std::to_string(n+1);
      inputs.set_value(pname, 1.0, info);
      if (info != 0) {
        uniform_projection_ = true;
        break;
      }
    }
  }

  // Set up the projector & add to variational parameters
  std::string strval, pname;
  if (uniform_projection_) {
    pjn_t pjn;
    gw_factor_.resize(1);

    // read projection type
    strval = inputs.set_value("gw_projection", "DOUBLON", info);
    boost::to_upper(strval);
    if (strval == "DOUBLON") {
      pjn = pjn_t::DOUBLON;
    }
    else if (strval == "HOLON") {
      pjn = pjn_t::HOLON;
    }
    else if (strval == "NONE") {
      pjn = pjn_t::NONE;
    }
    else {
      throw std::range_error("GW_Projector::init: invalid projection");
    }
    // set same projection at all sites
    for (auto& p : basis_projection_) p = pjn;
    for (auto& p : site_projection_) p = pjn;

    // read gw_factor
    if (pjn != pjn_t::NONE) {
      gw_factor_[0] = inputs.set_value("gw_factor", 1.0); 
    }
    else gw_factor_[0] = 1.0; 
    if (gw_factor_[0] < gw_cutoff()) {
      throw std::range_error("GW_Projector::init: `gw_factor' value out of range");
    }
    if (basis_projection_[0] != pjn_t::NONE) {
      vparms.add("gw_factor", gw_factor_[0], gw_cutoff(), 10.0);
    }
  }

  // site-type specific projection
  else {
    gw_factor_.resize(num_site_types_);

    for (int n=0; n<num_site_types_; ++n) {
      // read projection type
      pname = "gw_projection"+std::to_string(n+1);
      strval = inputs.set_value(pname, "DOUBLON", info);
      boost::to_upper(strval);
      if (strval == "DOUBLON") {
        basis_projection_[n] = pjn_t::DOUBLON;
      }
      else if (strval == "HOLON") {
        basis_projection_[n] = pjn_t::HOLON;
      }
      else if (strval == "NONE") {
        basis_projection_[n] = pjn_t::NONE;
      }
      else {
        throw std::range_error("GW_Projector::init: invalid projection");
      }

      // read gw_factor
      pname = "gw_factor"+std::to_string(n+1);
      if (basis_projection_[n] != pjn_t::NONE) {
        gw_factor_[n] = inputs.set_value(pname, 1.0);
      }
      else gw_factor_[n] = 1.0; 
      if (gw_factor_[n] < gw_cutoff()) {
        throw std::range_error("GW_Projector::init: 'gw_factor' value out of range");
      }
      if (basis_projection_[n] != pjn_t::NONE) {
        vparms.add(pname, gw_factor_[n], gw_cutoff(), 10.0);
      }
    }
    // set site-type specific projection 
    for (int i=0; i<num_sites_; ++i) {
      site_projection_[i] = basis_projection_[site_typeid_[i]];
    }
  }

  // 'default_case' true if  present & uniform & DOUBLON
  if (is_present_ && uniform_projection_ && site_projection_[0]==pjn_t::DOUBLON) {
    default_case_ = true;
  }
  else {
    default_case_ = false;
  }

  // GW ratio table
  set_ratio_table();

  return gw_factor_.size();
}

void GW_Projector::switch_off(void)
{
  is_present_ = false;
  default_case_ = false;
}

void GW_Projector::set_ratio_table(void)
{
  if (nickelate_2B_) {
    throw std::range_error("GW_Projector::set_ratio_table: case not implemented");
  }
  if (default_case_) {
    double g = gw_factor_[0];
    if ( g<0.0 ) {
      throw std::range_error("GW_Projector::set_ratio_table: out of range 'gw_factor' value");
    }
    gw_ratio_(0,0) = 1.0/(g*g); // nd_increament = -2
    gw_ratio_(0,1) = 1.0/g; // nd_increament = -1
    gw_ratio_(0,2) = 1.0; // nd_increament =  0
    gw_ratio_(0,3) = g; // nd_increament =  1
    gw_ratio_(0,4) = g*g; // nd_increament =  2
  }
  else if (uniform_projection_) {
    double g = gw_factor_[0];
    if ( g<0.0 ) {
      throw std::range_error("GW_Projector::set_ratio_table: out of range 'gw_factor' value");
    }
    // ratio
    for (int i=0; i<num_sites_; ++i) {
      gw_ratio_(i,0) = 1.0/(g*g); // nd_increament = -2
      gw_ratio_(i,1) = 1.0/g; // nd_increament = -1
      gw_ratio_(i,2) = 1.0; // nd_increament =  0
      gw_ratio_(i,3) = g; // nd_increament =  1
      gw_ratio_(i,4) = g*g; // nd_increament =  2
    }
  }
  else {
    // non-uniform projecton
    // sanity check
    for (int n=0; n<num_site_types_; ++n) {
      if (gw_factor_[n]<0.0 ){
        throw std::range_error("GW_Projector::set_ratio_table: out of range 'gw_factor' value");
      }
    }
    // set the ratio
    for (int i=0; i<num_sites_; ++i) {
      double g = gw_factor_[site_typeid_[i]];
      gw_ratio_(i,0) = 1.0/(g*g); // nd_increament = -2
      gw_ratio_(i,1) = 1.0/g; // nd_increament = -1
      gw_ratio_(i,2) = 1.0; // nd_increament =  0
      gw_ratio_(i,3) = g; // nd_increament =  1
      gw_ratio_(i,4) = g*g; // nd_increament =  2
    }
  }
}

bool GW_Projector::is_strong(void) const
{
  if (uniform_projection_) {
    return gw_factor_[0] <= 10.0*gw_cutoff();
  }
  else {
    for (int n=0; n<num_site_types_; ++n) {
      if (gw_factor_[n] <= 10.0*gw_cutoff()) return true;
    }
    return false;
  }
}

int GW_Projector::update_parameters(const VariationalParms& vparms)
{
  if (!is_present_) return 0;
  if (nickelate_2B_) {
    for (int n=0; n<3; ++n) {
      gw_factor_[n] = vparms[gfactor_names_[n]].value();
    }
    gfactor_ud_(0,0) = gw_factor_[0];
    gfactor_ud_(1,1) = gw_factor_[1];
    gfactor_ud_(0,1) = gw_factor_[2];
    gfactor_ud_(1,0) = gw_factor_[2];
    //gfactor_uu_(0,1) = gw_factor_[3];
    //gfactor_uu_(1,0) = gw_factor_[3];
    // ratio_table not used in this case
    return 0;
  }

  if (uniform_projection_) {
    if (basis_projection_[0] != pjn_t::NONE) {
      gw_factor_[0] = vparms["gw_factor"].value();
    }
  }
  else {
    for (int n=0; n<num_site_types_; ++n) {
      if (basis_projection_[n] != pjn_t::NONE) {
        gw_factor_[n] = vparms["gw_factor"+std::to_string(n+1)].value();
      }
    }
  }
  set_ratio_table();
  return 0;
}

double GW_Projector::gw_ratio(const vmc::BasisState& state, 
  const int& fr_site, const int& to_site, const vmc::move_t& move) const
{
  double gw_ratio = 1.0;
  if (!is_present_) return gw_ratio;
  if (fr_site == to_site) return gw_ratio;

  if (nickelate_2B_) {
    // upspin hop
    if (move == vmc::move_t::upspin_hop) {
      auto orbitals = state.parent_atom(fr_site);
      int m = site_typeid_[fr_site];
      for (const auto& s: orbitals) {
        int n = site_typeid_[s];
        if (state.op_ni_dn(s)) {
          gw_ratio /= gfactor_ud_(m,n);
          //std::cout << "g("<<m<<","<<n<<") = "<< gfactor_ud_(m,n) << "\n";
          //getchar();
        }
        if ((s!=fr_site) && state.op_ni_up(s)) {
          //gw_ratio /= gfactor_uu_(m,n); 
          gw_ratio /= gfactor_ud_(m,n); 
          //std::cout << "g("<<m<<","<<n<<") = "<< gfactor_ud_(m,n) << "\n";
          //getchar();
        }
      }

      // effect at the to_site
      orbitals = state.parent_atom(to_site);
      m = site_typeid_[to_site];
      for (const auto& s: orbitals) {
        int n = site_typeid_[s];
        if (state.op_ni_dn(s)) {
          gw_ratio *= gfactor_ud_(m,n);
        }
        if ((s!=fr_site) && state.op_ni_up(s)) {
          // if s==fr_site, the spin would have moved
          //gw_ratio *= gfactor_uu_(m,n); 
          gw_ratio *= gfactor_ud_(m,n); 
        }
        n++;
      }
    }

    // dnspin-hop
    else if (move == vmc::move_t::dnspin_hop) {
      auto orbitals = state.parent_atom(fr_site);
      int m = site_typeid_[fr_site];
      for (const auto& s: orbitals) {
        int n = site_typeid_[s];
        if (state.op_ni_up(s)) {
          gw_ratio /= gfactor_ud_(m,n);
        }
        if ((s!=fr_site) && state.op_ni_dn(s)) {
          //gw_ratio /= gfactor_uu_(m,n); 
          gw_ratio /= gfactor_ud_(m,n); 
        }
        n++;
      }

      // effect at the to_site
      orbitals = state.parent_atom(to_site);
      m = site_typeid_[to_site];
      for (const auto& s: orbitals) {
        int n = site_typeid_[s];
        if (state.op_ni_up(s)) {
          gw_ratio *= gfactor_ud_(m,n);
        }
        if ((s!=fr_site) && state.op_ni_dn(s)) {
          // if s==fr_site, the spin would have moved
          //gw_ratio *= gfactor_uu_(m,n); 
          gw_ratio *= gfactor_ud_(m,n); 
        }
        n++;
      }
    }
    else {
      throw std::range_error("GW_Projector::gw_ratio: undefined hopping move");
    }
    return gw_ratio;
  }


  if (default_case_) {
    if (state.op_ni_dblon(fr_site)) {
      // one doublon to be annihilated (change = -1)
      gw_ratio *= gw_ratio_(0,1);
    }
    if (!state.op_ni_holon(to_site)) {
      // one doublon to be created (change = +1)
      gw_ratio *= gw_ratio_(0,3);
    }
    return gw_ratio;
  }

  // fr_site
  pjn_t pjn = site_projection_[fr_site];
  if (pjn == pjn_t::DOUBLON) {
    if (state.op_ni_dblon(fr_site)) {
      // one doublon to be annihilated (change = -1)
      gw_ratio *= gw_ratio_(fr_site,1);
    }
  }
  else if (pjn == pjn_t::HOLON) {
    if (!state.op_ni_dblon(fr_site)) {
      // one holon to be created (change = 1)
      gw_ratio *= gw_ratio_(fr_site,3);
    }
  }

  // to_site
  pjn = site_projection_[to_site];
  if (pjn == pjn_t::DOUBLON) {
    if (!state.op_ni_holon(to_site)) {
      // one doublon to be created (change = +1)
      gw_ratio *= gw_ratio_(to_site,3);
    }
  }
  else if (pjn == pjn_t::HOLON) {
    if (state.op_ni_holon(to_site)) {
      // one holon to be annihilated (change = -1)
      gw_ratio *= gw_ratio_(to_site,1);
    }
  }
  return gw_ratio;

}

double GW_Projector::gw_ratio(const vmc::BasisState& state, 
  const int& fr_site, const int& to_site) const
{
  double gw_ratio = 1.0;
  if (!is_present_) return gw_ratio;
  if (fr_site == to_site) return gw_ratio;

  if (default_case_) {
    if (state.op_ni_dblon(fr_site)) {
      // one doublon to be annihilated (change = -1)
      gw_ratio *= gw_ratio_(0,1);
    }
    if (!state.op_ni_holon(to_site)) {
      // one doublon to be created (change = +1)
      gw_ratio *= gw_ratio_(0,3);
    }
    return gw_ratio;
  }

  // fr_site
  pjn_t pjn = site_projection_[fr_site];
  if (pjn == pjn_t::DOUBLON) {
    if (state.op_ni_dblon(fr_site)) {
      // one doublon to be annihilated (change = -1)
      gw_ratio *= gw_ratio_(fr_site,1);
    }
  }
  else if (pjn == pjn_t::HOLON) {
    if (!state.op_ni_dblon(fr_site)) {
      // one holon to be created (change = 1)
      gw_ratio *= gw_ratio_(fr_site,3);
    }
  }

  // to_site
  pjn = site_projection_[to_site];
  if (pjn == pjn_t::DOUBLON) {
    if (!state.op_ni_holon(to_site)) {
      // one doublon to be created (change = +1)
      gw_ratio *= gw_ratio_(to_site,3);
    }
  }
  else if (pjn == pjn_t::HOLON) {
    if (state.op_ni_holon(to_site)) {
      // one holon to be annihilated (change = -1)
      gw_ratio *= gw_ratio_(to_site,1);
    }
  }
  return gw_ratio;
}

double GW_Projector::gw_ratio_pairhop(const int& fr_site, const int& to_site) const
{
  double gw_ratio = 1.0;
  if (default_case_) return gw_ratio;
  if (!is_present_) return gw_ratio;
  if (uniform_projection_) return gw_ratio;
  if (fr_site == to_site) return gw_ratio;
  if (nickelate_2B_) {
    throw std::range_error("GW_Projector::gw_ratio_pairhop: case not implemented");
  }

  // for non-uniform projection
  // fr_site
  pjn_t pjn = site_projection_[fr_site];
  if (pjn == pjn_t::DOUBLON) {
    // one doublon to be annihilated (change = -1)
    gw_ratio *= gw_ratio_(fr_site,1);
  }
  else if (pjn == pjn_t::HOLON) {
    // one holon to be created (change = 1)
    gw_ratio *= gw_ratio_(fr_site,3);
  }

  // to_site
  pjn = site_projection_[to_site];
  if (pjn == pjn_t::DOUBLON) {
    // one doublon to be created (change = +1)
    gw_ratio *= gw_ratio_(to_site,3);
  }
  else if (pjn == pjn_t::HOLON) {
    // one holon to be annihilated (change = -1)
    gw_ratio *= gw_ratio_(to_site,1);
  }
  return gw_ratio;
}

void GW_Projector::get_grad_logp(const vmc::BasisState& state, RealVector& grad) const
{
  if (!is_present_) return;
  if (default_case_) {
    int nd = 0;
    for (int n=0; n<num_sites_; ++n) {
      nd += state.op_ni_dblon(n);
    } 
    grad(0) = static_cast<double>(nd)/gw_factor_[0];
    return;
  }

  if (nickelate_2B_) {
    // total no of DOUBLON in each orbital
    int nd_ss = 0;
    int nd_dd = 0;
    int nd_sd = 0;
    //int nd_uu = 0;
    for (int i=0; i<state.num_atoms(); ++i) {
      auto orbitals = state.atom(i);
      for (const auto& s : orbitals) {
        int n = state.op_ni_dblon(s);
        if (site_typeid_[s]==0) nd_ss += n;
        if (site_typeid_[s]==1) nd_dd += n;
      }
      if (state.op_ni_up(orbitals[0]) && state.op_ni_dn(orbitals[1])) {
        nd_sd += 1;  
      }
      if (state.op_ni_dn(orbitals[0]) && state.op_ni_up(orbitals[1])) {
        nd_sd += 1;  
      }
      if (state.op_ni_up(orbitals[0]) && state.op_ni_up(orbitals[1])) {
        //nd_uu += 1;  
        nd_sd += 1;  
      }
      if (state.op_ni_dn(orbitals[0]) && state.op_ni_dn(orbitals[1])) {
        //nd_uu += 1;  
        nd_sd += 1;  
      }
    }
    grad(0) = static_cast<double>(nd_ss)/gw_factor_[0];
    grad(1) = static_cast<double>(nd_dd)/gw_factor_[1];
    grad(2) = static_cast<double>(nd_sd)/gw_factor_[2];
    //grad(3) = static_cast<double>(nd_uu)/gw_factor_[3];

    return;
  }



  if (uniform_projection_) {
    int nd = 0;
    pjn_t pjn = site_projection_[0]; 
    if (pjn == pjn_t::DOUBLON) {
      for (int n=0; n<num_sites_; ++n) {
        nd += state.op_ni_dblon(n);
      } 
    }
    else if (pjn == pjn_t::HOLON) {
      for (int n=0; n<num_sites_; ++n) {
        nd += state.op_ni_holon(n);
      }  
    }
    grad(0) = static_cast<double>(nd)/gw_factor_[0];
  }
  else {
    // For each sublattice: total no of DOUBLON or HOLON
    std::vector<int> N_total(num_site_types_,0);
    for (int i=0; i<num_sites_; ++i) {
      pjn_t pjn = site_projection_[i];
      int stype = site_typeid_[i];
      if (pjn == pjn_t::DOUBLON) {
        N_total[stype] += state.op_ni_dblon(i);
      }
      else if (pjn == pjn_t::HOLON) {
        N_total[stype] += state.op_ni_holon(i);
      }
    }
    // Gradient. (in case of pjn_t::NONE, N_total[n]=0, hence grad[n]=0)
    for (int n=0; n<num_site_types_; ++n) {
      grad(n) = static_cast<double>(N_total[n])/gw_factor_[n];
    }
  }
}

/*--------------------------------------------------
* WavefunProjector: Handles all projectors
*---------------------------------------------------*/
void WavefunProjector::init(const lattice::Lattice& lattice, const input::Parameters& inputs) 
{
  num_sites_ = lattice.num_sites();

  varparms_.clear();

  // GW projector
  gw_projector_.init(lattice, inputs, varparms_);

  // DH projector
  init_dh_projector(lattice, inputs);
}


void WavefunProjector::update(const input::Parameters& inputs) 
{ 
  for (auto& p : varparms_) {
    double x = inputs.set_value(p.name(), p.value());
    p.change_value(x);
  }

  // update GW parameters
  gw_projector_.update_parameters(varparms_);

  // dh factors
  if (dh_projector_) {
    if (dh_range_ >= 1) dh_factor1_ = varparms_["dhfactor1"].value();
    if (dh_range_ >= 2) dh_factor2_ = varparms_["dhfactor2"].value();
    if (dh_range_ >= 3) dh_factor3_ = varparms_["dhfactor3"].value();
  }
} 

void WavefunProjector::update(const var::parm_vector& pvector, const unsigned& start_pos)
{ 
  //varparms_.update(vparms,begin,end); 
  unsigned i = 0;
  for (auto& p : varparms_) {
    p.change_value(pvector[start_pos+i]);
    i++;
  }

  // update GW parameters
  gw_projector_.update_parameters(varparms_);

  // dh factors
  if (dh_projector_) {
    if (dh_range_ >= 1) dh_factor1_ = varparms_["dhfactor1"].value();
    if (dh_range_ >= 2) dh_factor2_ = varparms_["dhfactor2"].value();
    if (dh_range_ >= 3) dh_factor3_ = varparms_["dhfactor3"].value();
  }
}

void WavefunProjector::get_vparm_names(std::vector<std::string>& vparm_names, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_names[start_pos+i] = p.name(); ++i;
  }
}

void WavefunProjector::get_vparm_values(var::parm_vector& vparm_values, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void WavefunProjector::get_vparm_vector(std::vector<double>& vparm_values, unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_values[start_pos+i] = p.value(); ++i;
  }
}

void WavefunProjector::get_vparm_lbound(var::parm_vector& vparm_lb, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_lb[start_pos+i] = p.lbound(); ++i;
  }
}

void WavefunProjector::get_vparm_ubound(var::parm_vector& vparm_ub, 
  unsigned start_pos) const
{
  unsigned i = 0;
  for (auto& p : varparms_) {
    vparm_ub[start_pos+i] = p.ubound(); ++i;
  }
}

double WavefunProjector::dh_factor1(void) const 
{
  return varparms_["dh_factor1"].value();
}

double WavefunProjector::dh_factor2(void) const 
{
  return varparms_["dh_factor2"].value();
}

double WavefunProjector::dh_factor3(void) const 
{
  return varparms_["dh_factor3"].value();
}


int WavefunProjector::init_dh_projector(const lattice::Lattice& lattice, 
  const input::Parameters& inputs) 
{
  int info;
  dh_projector_ = inputs.set_value("dh_projector",false, info);
  if (!dh_projector_) return 0;

  // dh parameters
  dh_range_ = 1;
  // nn-1
  dh_factor1_ = inputs.set_value("dhfactor1", 1.0);
  if (dh_factor1_<0.0) throw std::range_error("WavefunProjector::init: out-of-range 'dh_factor'-value");
  varparms_.add("dhfactor1", dh_factor1_, 1.0E-4, 2.0);
  // nn-2
  dh_factor2_ = inputs.set_value("dhfactor2", 1.0, info);
  if (info == 0) {
    if (dh_factor2_<0.0) throw std::range_error("WavefunProjector::init: out-of-range 'dh_factor'-value");
    varparms_.add("dhfactor2", dh_factor2_, 1.0E-4, 2.0);
    dh_range_++;
    // nn-3
    dh_factor3_ = inputs.set_value("dhfactor3", 1.0, info);
    if (info == 0) {
      if (dh_factor3_<0.0) throw std::range_error("WavefunProjector::init: out-of-range 'dh_factor'-value");
      varparms_.add("dhfactor3", dh_factor3_, 1.0E-4, 2.0);
      dh_range_++;
    }
  }

  // Create Adjacency table

  // determine 3 smallest NN distances in the lattice
  std::vector<double> dist;
  for (int i=1; i<num_sites_; ++i) {
    dist.push_back(lattice.site(i).coord().norm());
  }
  std::sort(dist.begin(),dist.end());
  std::vector<double> nn_dist;
  nn_dist.push_back(dist[0]);
  for (const auto& d : dist) {
    if (d > nn_dist.back()) nn_dist.push_back(d);
    if (nn_dist.size() >= 3) break;
  }
  /*
  std::cout << "dnn1 = " << nn_dist[0] << "\n";
  std::cout << "dnn2 = " << nn_dist[1] << "\n";
  std::cout << "dnn3 = " << nn_dist[2] << "\n";
  */

  // lattice size
  std::vector<double> Lsize(3);
  double L1 = lattice.size1()*lattice.vector_a1().norm();
  double L2 = lattice.size2()*lattice.vector_a2().norm();
  double L3 = lattice.size3()*lattice.vector_a3().norm();

  // Adjacency lists
  adjacency_.clear();
  adjacency_.resize(num_sites_);
  for (int i=0; i<num_sites_; ++i) {
    Vector3d Ri = lattice.site(i).coord();
    for (int j=0; j<num_sites_; ++j) {
      if (i == j) continue;
      Vector3d Rj = lattice.site(j).coord();
      Vector3d R = Ri-Rj;
      R = R.cwiseAbs();
      if (lattice.bc1() != lattice::boundary_type::open) {
        if (R[0] > L1/2) R[0] = L1-R[0];
      }
      if (lattice.bc2() != lattice::boundary_type::open) {
        if (R[1] > L2/2) R[1] = L2-R[1];
      }
      if (lattice.bc3() != lattice::boundary_type::open) {
        if (R[2] > L3/2) R[2] = L3-R[2];
      }
      double d = R.norm();
      //std::cout << i << "  " << j << "  " << d << "\n"; getchar();
      if (std::abs(d-nn_dist[0])<1.0E-6) adjacency_[i].add_nn1site(j);
      if (std::abs(d-nn_dist[1])<1.0E-6) adjacency_[i].add_nn2site(j);
      if (std::abs(d-nn_dist[2])<1.0E-6) adjacency_[i].add_nn3site(j);
    }
  }
  // check 
  /*
  for (int i=0; i<num_sites_; ++i) {
    std::cout << i << ": ";
    for (const auto& nn : adjacency_[i].nn1sites() ) {
      std::cout << nn << ", ";
    }
    std::cout << "\n";
  }
  getchar();
  */
  return 0;
}

double WavefunProjector::dh_ratio(const vmc::BasisState& state, 
  const int& fr_site, const int& to_site) const
{
  if (!dh_projector_) return 1.0;
  if (fr_site == to_site) return 1.0;

  double dh_ratio = 1.0;
  int nholon, ndblon;
  // If 'fr_site' site was a doublon
  if (state.op_ni_dblon(fr_site)) {
    // Effect at 1st NN
    // if it was isolated in 1st NN
    nholon = 0;
    for (const auto& nn : adjacency_[fr_site].nn1sites()) {
      if (state.op_ni_holon(nn)) {
        nholon++;
        // if it makes this holon isolated
        ndblon = 0; 
        for (const auto& mm : adjacency_[nn].nn1sites()) {
          ndblon += state.op_ni_dblon(mm);
        }
        if (ndblon == 1) dh_ratio *= dh_factor1_; // isolated holon created
      }
    }
    if (nholon == 0) dh_ratio /= dh_factor1_; // isolated doublon removed

    // Effect at 2nd NN
    if (dh_range_>=2) {
      nholon = 0;
      for (const auto& nn : adjacency_[fr_site].nn2sites()) {
        if (state.op_ni_holon(nn)) {
          nholon++;
          // if it makes this holon isolated
          ndblon = 0; 
          for (const auto& mm : adjacency_[nn].nn2sites()) {
            ndblon += state.op_ni_dblon(mm);
          }
          if (ndblon == 1) dh_ratio *= dh_factor2_; // isolated holon created
        }
      }
      if (nholon == 0) dh_ratio /= dh_factor2_; // isolated doublon removed
    }

    // Effect at 3rd NN
    if (dh_range_>=3) {
      nholon = 0;
      for (const auto& nn : adjacency_[fr_site].nn3sites()) {
        if (state.op_ni_holon(nn)) {
          nholon++;
          // if it makes this holon isolated
          ndblon = 0; 
          for (const auto& mm : adjacency_[nn].nn3sites()) {
            ndblon += state.op_ni_dblon(mm);
          }
          if (ndblon == 1) dh_ratio *= dh_factor3_; // isolated holon created
        }
      }
      if (nholon == 0) dh_ratio /= dh_factor3_; // isolated doublon removed
    }
  }
  else {
    // A holon is created at the 'fr_site'
    // Effect at 1st NN - if it is isolated
    ndblon = 0;
    for (const auto& nn : adjacency_[fr_site].nn1sites()) {
      if (state.op_ni_dblon(nn)) {
        ndblon++;
        // if it was an isolated doublon
        nholon = 0; 
        for (const auto& mm : adjacency_[nn].nn1sites()) {
          nholon += state.op_ni_holon(mm);
        }
        if (nholon == 0) dh_ratio /= dh_factor1_; // doublon isolation removed
      }
    }
    if (ndblon == 0) dh_ratio *= dh_factor1_; // isolated holon created

    // Effect at 2nd NN 
    if (dh_range_>=2) {
      ndblon = 0;
      for (const auto& nn : adjacency_[fr_site].nn2sites()) {
        if (state.op_ni_dblon(nn)) {
          ndblon++;
          // if it was an isolated doublon
          nholon = 0; 
          for (const auto& mm : adjacency_[nn].nn2sites()) {
            nholon += state.op_ni_holon(mm);
          }
          if (nholon == 0) dh_ratio /= dh_factor2_; // doublon isolation removed
        }
      }
      if (ndblon == 0) dh_ratio *= dh_factor2_; // isolated holon created
    }

    // Effect at 3rd NN 
    if (dh_range_>=3) {
      ndblon = 0;
      for (const auto& nn : adjacency_[fr_site].nn3sites()) {
        if (state.op_ni_dblon(nn)) {
          ndblon++;
          // if it was an isolated doublon
          nholon = 0; 
          for (const auto& mm : adjacency_[nn].nn3sites()) {
            nholon += state.op_ni_holon(mm);
          }
          if (nholon == 0) dh_ratio /= dh_factor3_; // doublon isolation removed
        }
      }
      if (ndblon == 0) dh_ratio *= dh_factor3_; // isolated holon created
    }
  }

  // If 'to_site' was a holon
  bool new_holon;
  if (state.op_ni_holon(to_site)) {
    // Effect at 1st NN
    // if it was isolated in 1st NN
    ndblon = 0;
    for (const auto& nn : adjacency_[to_site].nn1sites()) {
      if (nn==fr_site) continue; // even if it was doublon, it was annilated
      if (state.op_ni_dblon(nn)) {
        ndblon++;
        // if it makes this doublon isolated
        nholon = 0; 
        for (const auto& mm : adjacency_[nn].nn1sites()) {
          if (mm==fr_site) {
            if (!state.op_ni_dblon(mm)) nholon++; // a holon was created here
          } 
          else {
            nholon += state.op_ni_holon(mm);
          }
        }
        if (nholon == 1) dh_ratio *= dh_factor1_; // doublon got isolated
      }
    }
    if (ndblon == 0) dh_ratio /= dh_factor1_; // isolated holon removed

    // Effect at 2nd NN
    // if it was isolated in 2nd NN
    if (dh_range_>=2) {
      ndblon = 0;
      for (const auto& nn : adjacency_[to_site].nn2sites()) {
        if (nn==fr_site) continue; // even if it was doublon, it was annilated
        if (state.op_ni_dblon(nn)) {
          ndblon++;
          // if it makes this doublon isolated
          nholon = 0; 
          for (const auto& mm : adjacency_[nn].nn2sites()) {
            if (mm==fr_site) {
              if (!state.op_ni_dblon(mm)) nholon++; // a holon was created here
            } 
            else {
              nholon += state.op_ni_holon(mm);
            }
          }
          if (nholon == 1) dh_ratio *= dh_factor2_; // isolated doublon created
        }
      }
      if (ndblon == 0) dh_ratio /= dh_factor2_; // isolated holon removed
    }

    // Effect at 3rd NN
    // if it was isolated in 3rd NN
    if (dh_range_>=3) {
      ndblon = 0;
      for (const auto& nn : adjacency_[to_site].nn3sites()) {
        if (nn==fr_site) continue; // even if it was doublon, it was annilated
        if (state.op_ni_dblon(nn)) {
          ndblon++;
          // if it makes this doublon isolated
          nholon = 0; 
          for (const auto& mm : adjacency_[nn].nn3sites()) {
            if (mm==fr_site) {
              if (!state.op_ni_dblon(mm)) nholon++; // a holon was created here
            } 
            else {
              nholon += state.op_ni_holon(mm);
            }
          }
          if (nholon == 1) dh_ratio *= dh_factor3_; // isolated doublon created
        }
      }
      if (ndblon == 0) dh_ratio /= dh_factor3_; // isolated holon removed
    }
  }
  else {
    // A doublon is created at the 'to_site'
    // Effect at 1st NN - if it is isolated
    nholon = 0;
    for (const auto& nn : adjacency_[to_site].nn1sites()) {
      if (nn==fr_site && !state.op_ni_dblon(nn)) new_holon = true;
      else new_holon = false;
      if (state.op_ni_holon(nn) || new_holon) {
        nholon++;
        // if it was an isolated holon
        ndblon = 0; 
        for (const auto& mm : adjacency_[nn].nn1sites()) {
          if (mm==fr_site) continue; // even if doublon, it was removed
          ndblon += state.op_ni_dblon(mm);
        }
        if (ndblon == 0) dh_ratio /= dh_factor1_; // doublon isolation removed
      }
    }
    if (nholon == 0) dh_ratio *= dh_factor1_; // isolated doublon created

    if (dh_range_>=2) {
      nholon = 0;
      for (const auto& nn : adjacency_[to_site].nn2sites()) {
        if (nn==fr_site && !state.op_ni_dblon(nn)) new_holon = true;
        else new_holon = false;
        if (state.op_ni_holon(nn) || new_holon) {
          nholon++;
          // if it was an isolated holon
          ndblon = 0; 
          for (const auto& mm : adjacency_[nn].nn2sites()) {
            if (mm==fr_site) continue; // even if doublon, it was removed
            ndblon += state.op_ni_dblon(mm);
          }
          if (ndblon == 0) dh_ratio /= dh_factor2_; // doublon isolation removed
        }
      }
      if (nholon == 0) dh_ratio *= dh_factor2_; // isolated doublon created
    }

    if (dh_range_>=3) {
      nholon = 0;
      for (const auto& nn : adjacency_[to_site].nn3sites()) {
        if (nn==fr_site && !state.op_ni_dblon(nn)) new_holon = true;
        else new_holon = false;
        if (state.op_ni_holon(nn) || new_holon) {
          nholon++;
          // if it was an isolated holon
          ndblon = 0; 
          for (const auto& mm : adjacency_[nn].nn3sites()) {
            if (mm==fr_site) continue; // even if doublon, it was removed
            ndblon += state.op_ni_dblon(mm);
          }
          if (ndblon == 0) dh_ratio /= dh_factor3_; // doublon isolation removed
        }
      }
      if (nholon == 0) dh_ratio *= dh_factor3_; // isolated doublon created
    }
  }

  return dh_ratio;
}


void WavefunProjector::get_grad_logp(const vmc::BasisState& state, RealVector& grad) const
{
  gw_projector_.get_grad_logp(state, grad);
  /*
  // Calculate Gradient-of-log-of-projectors
  // GW projector
  if (gw_factor_.size() == 1) {
    // uniform projection case
    int num_dblon = 0;
    for (int n=0; n<num_sites_; ++n) {
      if (state.op_ni_dblon(n)) num_dblon += 1;
    } 
    grad(0) = static_cast<double>(num_dblon)/gw_factor_[0];
  }
  else {
    // sublattice-selective projection case
    std::vector<int> num_dblon(num_site_types_,0);
    std::vector<int> num_holon(num_site_types_,0);
    for (int n=0; n<num_sites_; ++n) {
      int sublatt = sitetype_list_[n];
      if (gw_ptypes_[sublatt] == gw_ptype::DOUBLON) {
        if (state.op_ni_dblon(n)) num_dblon[sublatt] += 1;
      }
      if (gw_ptypes_[sublatt] == gw_ptype::HOLON) {
        if (state.op_ni_holon(n)) num_holon[sublatt] += 1;
      }
    }
    for (int n=0; n<num_site_types_; ++n) {
      if (gw_ptypes_[n] == gw_ptype::DOUBLON) {
        grad(n) = static_cast<double>(num_dblon[n])/gw_factor_[n];
      }
      else if (gw_ptypes_[n] == gw_ptype::HOLON) {
        grad(n) = static_cast<double>(num_holon[n])/gw_factor_[n];
      }
      else {
        throw std::logic_error("WavefunProjector::get_grad_logp: logic error");
      }
    }
  }
  */

  // DH projector
  if (dh_projector_) {
    throw std::range_error("WavefunProjector::get_grad_logp: not implemented for DH projector\n");
  }
}


} // end namespace var

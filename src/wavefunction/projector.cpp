/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-16 23:17:49
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-03-05 11:47:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./projector.h"

namespace var {

void WavefunProjector::init(const lattice::Lattice& lattice, const input::Parameters& inputs) 
{
  num_sites_ = lattice.num_sites();

  /*
  int info;
  gw_projector_ = inputs.set_value("gw_projector", false, info);
  double g = 1.0;
  if (gw_projector_) {
    g = inputs.set_value("gfactor", 1.0);
    if (g<0.0) throw std::range_error("WavefunProjector::init: out-of-range 'g'-value.");
    varparms_.add("gfactor", g, 1.0E-4, 10.0);
  }
  //pfactors_.insert({"g", g});
  // gw ratio
  gw_ratio_.resize(5);
  gw_ratio_[0] = 1.0; // nd_increament = -2
  gw_ratio_[1] = 1.0; // nd_increament = -1
  gw_ratio_[2] = 1.0; // nd_increament =  0
  gw_ratio_[3] = 1.0; // nd_increament =  1
  gw_ratio_[4] = 1.0; // nd_increament =  2
  */

  // GW projector
  init_gw_projector(lattice, inputs);

  // DH projector
  init_dh_projector(lattice, inputs);
}

void WavefunProjector::update(const input::Parameters& inputs) 
{ 
  for (auto& p : varparms_) {
    double x = inputs.set_value(p.name(), p.value());
    p.change_value(x);
  }

  if (gw_projector_) {
    if (num_site_types_==1) {
      gw_factor_[0] = varparms_["gw_factor"].value();
    }
    else {
      for (int n=0; n<num_site_types_; ++n) {
        gw_factor_[n] = varparms_["gw_factor"+std::to_string(n+1)].value();
      }
    }
    set_gw_ratio();
  }

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

  // gw factors
  if (gw_projector_) {
    if (num_site_types_==1) {
      gw_factor_[0] = varparms_["gw_factor"].value();
    }
    else {
      for (int n=0; n<num_site_types_; ++n) {
        gw_factor_[n] = varparms_["gw_factor"+std::to_string(n+1)].value();
      }
    }
    set_gw_ratio();
  }

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

bool WavefunProjector::gwfactor_is_zero(void) const
{
  if (num_site_types_ == 1) {
    return varparms_["gw_factor"].value() <= gw_cutoff();
  }
  else {
    return false;
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

void WavefunProjector::get_grad_logp(const vmc::BasisState& state, RealVector& grad) const
{
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
    for (int n=0; n<num_sites_; ++n) {
      int type = sitetype_list_[n];
      if (state.op_ni_dblon(n)) num_dblon[type] += 1;
    }
    for (int i=0; i<num_site_types_; ++i) {
      grad(i) = static_cast<double>(num_dblon[i])/gw_factor_[i];
    }
  }

  // DH projector
  if (dh_projector_) {
    throw std::range_error("WavefunProjector::get_grad_logp: not implemented for DH projector\n");
  }
}



int WavefunProjector::init_gw_projector(const lattice::Lattice& lattice, const input::Parameters& inputs) 
{

  num_site_types_ = lattice.num_site_types();
  sitetype_list_.resize(num_sites_);
  for (const auto& s : lattice.sites()) sitetype_list_[s.id()] = s.type();

  if (lattice.id()==lattice::lattice_id::SQUARE_2SITE) {
    // different 'gfactor' at different sublattice
    num_site_types_ = lattice.num_site_types();
    gw_factor_.resize(num_site_types_);
    sitetype_list_.resize(num_sites_);
    for (const auto& s : lattice.sites()) sitetype_list_[s.id()] = s.type();
  }
  else {
    // uniform 'gfactor' 
    num_site_types_ = 1;
    gw_factor_.resize(num_site_types_);
    for (int i=0; i<num_sites_; ++i) sitetype_list_[i] = 0;
  }

  // read from input
  int info;
  gw_projector_ = inputs.set_value("gw_projector",false,info);
  if (gw_projector_) {
    for (int n=0; n<num_site_types_; ++n) {
      std::string pname;
      if (num_site_types_ == 1) pname = "gw_factor";
      else pname = "gw_factor"+std::to_string(n+1);
      double g = inputs.set_value(pname, 1.0);
      //if (g<0.0) throw std::range_error("WavefunProjector::init: out-of-range 'g'-value.");
      if (g<gw_cutoff()) throw std::range_error("WavefunProjector::init: out-of-range 'g'-value.");
      gw_factor_[n] = g;
      //varparms_.add(pname, g, 1.0E-4, 10.0);
      varparms_.add(pname, g, gw_cutoff(), 10.0);
    }
  }
  else {
    for (int n=0; n<num_site_types_; ++n) gw_factor_[n] = 1.0;
  }

  // set ratio table
  gw_ratio_table_.resize(num_site_types_);
  for (auto& elem : gw_ratio_table_) elem.resize(5);
  set_gw_ratio();

  return 0;
}

void WavefunProjector::set_gw_ratio(void) 
{ 
  for (int n=0; n<num_site_types_; ++n) {
    double g = gw_factor_[n];
    if (g<0.0) throw std::range_error("WavefunProjector::set_gw_ratio: out-of-range 'g'-factor");
    gw_ratio_table_[n][0] = 1.0/(g*g); // nd_increament = -2
    gw_ratio_table_[n][1] = 1.0/g; // nd_increament = -1
    gw_ratio_table_[n][2] = 1.0; // nd_increament =  0
    gw_ratio_table_[n][3] = g; // nd_increament =  1
    gw_ratio_table_[n][4] = g*g; // nd_increament =  2
  }
  /*
  double g = gw_factor();
  //if (g<0.0 || g >1.0) throw std::range_error("WavefunProjector::set_gw_ratio: out-of-range 'g'-value");
  if (g<0.0) throw std::range_error("WavefunProjector::set_gw_ratio: out-of-range 'g'-value");
  gw_ratio_[0] = 1.0/(g*g);  // nd_increament = -2
  gw_ratio_[1] = 1.0/g;      // nd_increament = -1
  gw_ratio_[2] = 1.0;        // nd_increament =  0
  gw_ratio_[3] = g;          // nd_increament =  1
  gw_ratio_[4] = g*g;        // nd_increament =  2
  */
} 


double WavefunProjector::gw_ratio(const vmc::BasisState& state, 
  const int& fr_site, const int& to_site) const
{
  if (!gw_projector_) return 1.0;
  if (fr_site == to_site) return 1.0;

  double gw_ratio = 1.0;
  if (state.op_ni_dblon(fr_site)) {
    // one doublon to be annihilated (change = -1)
    int stype = sitetype_list_[fr_site];
    gw_ratio *= gw_ratio_table_[stype][1];
  }

  if (!state.op_ni_holon(to_site)) {
    // one doublon to be created (change = +1)
    int stype = sitetype_list_[to_site];
    gw_ratio *= gw_ratio_table_[stype][3];
  }

  return gw_ratio;
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


} // end namespace var

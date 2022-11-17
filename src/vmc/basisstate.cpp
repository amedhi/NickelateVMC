/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:20:28
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-02-20 04:54:18
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "basisstate.h"
#include <stdexcept>
#include <algorithm>

namespace vmc {

void BasisState::init(const lattice::Lattice& lattice, const model::Hamiltonian& model)
{
  num_sites_ = lattice.num_sites();
  site_states_.resize(num_sites_);
  for (auto& state : site_states_) state.reset();

  // set model projection
  site_projection_.resize(num_sites_);
  if (model.projection()) {
    have_projection_ = true;
    for (int i=0; i<num_sites_; ++i) {
      int stype = lattice.site(i).type();
      switch (model.projection().get(stype)) {
        case model::projection_t::DOUBLON:
          site_projection_[i] = pjn_t::DOUBLON;
          break;
        case model::projection_t::HOLON:
          site_projection_[i] = pjn_t::HOLON;
          break;
        default:
          site_projection_[i] = pjn_t::NONE;
      }
    }
  }
  else {
    have_projection_ = false;
    for (int i=0; i<num_sites_; ++i) {
      site_projection_[i] = pjn_t::NONE;
    }
  }

  // clear
  upspin_sites_.clear();
  dnspin_sites_.clear();
  uphole_sites_.clear();
  dnhole_sites_.clear();

  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
}

void BasisState::clear(void)
{
  for (auto& state : site_states_) state.reset();
  upspin_sites_.clear();
  dnspin_sites_.clear();
  uphole_sites_.clear();
  dnhole_sites_.clear();
} 

void BasisState::init_spins(const int& num_upspins, const int& num_dnspins)
{
  // constraint checks
  if (num_upspins>num_sites_ || num_dnspins>num_sites_) {
    throw std::range_error(" BasisState::init_spins: spin number exceed capacity");
  }

  // max states considering DOUBLON projection
  int max_states = 0;
  for (const auto& p : site_projection_) {
    if (p == pjn_t::DOUBLON) max_states += 1;
    else max_states += 2;
  }
  if (num_upspins+num_dnspins>max_states) {
    throw std::range_error(" BasisState::init_spins: spin number exceed capacity");
  }

  // min particle number considering HOLON projection
  int min_particles = 0;
  for (const auto& p : site_projection_) {
    if (p == pjn_t::HOLON) min_particles += 1;
  }
  if (num_upspins+num_dnspins<min_particles) {
    throw std::range_error(" BasisState::init_spins: spin number not enough to meet projection");
  }

  // things fine, go ahead
  num_upspins_ = num_upspins;
  num_dnspins_ = num_dnspins;
  num_upholes_ = num_sites_ - num_upspins;
  num_dnholes_ = num_sites_ - num_dnspins;

  upspin_sites_.resize(num_upspins_);
  dnspin_sites_.resize(num_dnspins_);
  uphole_sites_.resize(num_upholes_);
  dnhole_sites_.resize(num_dnholes_);

  // random generators
  // rng site generator
  rng_.set_site_generator(0,num_sites_-1);
  // assuming 'num_upspins>0', 'num_dnspins>0'
  rng_.set_upspin_generator(0,num_upspins_-1);
  rng_.set_dnspin_generator(0,num_dnspins_-1);
  // hole numbers may be zero
  int m = std::max(num_upholes_,1);
  rng_.set_uphole_generator(0,m-1);
  int n = std::max(num_dnholes_,1);
  rng_.set_dnhole_generator(0,n-1);
} 

void BasisState::set_random(void)
{
  //set_random_old(); return;

  // reset state
  for (auto& state : site_states_)  state.reset();

  // HOLON projected sites must be filled first
  int next_up = 0;
  int next_dn = 0;
  std::vector<int> no_holon_sites;
  for (int n=0; n<num_sites_; ++n) {
    if (site_projection_[n]==pjn_t::HOLON) {
      no_holon_sites.push_back(n);
    }
  }
  if (no_holon_sites.size() > 0) {
    std::shuffle(no_holon_sites.begin(),no_holon_sites.end(),rng_);
    // UP spins & holes
    int m;
    int idx = 0;
    for (m=0; m<num_upspins_; ++m) {
      if (idx >= no_holon_sites.size()) break;
      int site = no_holon_sites[idx];
      site_states_[site].put_upspin(m);
      upspin_sites_[m] = site;
      idx++;
    }
    next_up = m;
    // if all 'no_holon_sites' are not filled
    for (m=0; m<num_dnspins_; ++m) {
      if (idx >= no_holon_sites.size()) break;
      int site = no_holon_sites[idx];
      site_states_[site].put_dnspin(m);
      dnspin_sites_[m] = site;
      idx++;
    }
    next_dn = m;
  }

  // fill all other states
  std::vector<int> all_sites(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_sites[i] = i;
  std::shuffle(all_sites.begin(),all_sites.end(),rng_);

  // fill the remaining UP spins  
  for (const auto& site : all_sites) {
    if (next_up >= num_upspins_) break;
    if (site_states_[site].have_upspin()) continue;
    if (site_projection_[site]==pjn_t::DOUBLON 
     && site_states_[site].is_not_empty()) continue;
    site_states_[site].put_upspin(next_up);
    upspin_sites_[next_up] = site;
    next_up++;
  }

  // fill the remaining DN spins  
  std::shuffle(all_sites.begin(),all_sites.end(),rng_);
  for (const auto& site : all_sites) {
    if (next_dn >= num_dnspins_) break;
    if (site_states_[site].have_dnspin()) continue;
    if (site_projection_[site]==pjn_t::DOUBLON 
     && site_states_[site].is_not_empty()) continue;
    site_states_[site].put_dnspin(next_dn);
    dnspin_sites_[next_dn] = site;
    next_dn++;
  }

  // sanity check
  if (next_up!=num_upspins_ && next_dn!=num_dnspins_) {
    std::cout << next_up << "  " << num_upspins_ << "\n";
    std::cout << next_dn << "  " << num_dnspins_ << "\n";
    throw std::logic_error(" BasisState::set_random: ivalid state status-1");
  }
  int n_up=0;
  int n_dn=0;
  for (const auto& state : site_states_) {
    if (state.have_upspin()) n_up++;
    if (state.have_dnspin()) n_dn++;
  }
  if (n_up!=num_upspins_ && n_dn!=num_dnspins_) {
    throw std::logic_error(" BasisState::set_random: ivalid state status-2");
  }

  // allotment of 'holes' 
  next_up = 0;
  for (const auto& site : all_sites) {
    if (next_up >= num_upholes_) break;
    if (!site_states_[site].have_upspin()) {
      site_states_[site].put_uphole(next_up);
      uphole_sites_[next_up] = site;
      next_up++;
    }
  }
  next_dn = 0;
  for (const auto& site : all_sites) {
    if (next_dn >= num_dnholes_) break;
    if (!site_states_[site].have_dnspin()) {
      site_states_[site].put_dnhole(next_dn);
      dnhole_sites_[next_dn] = site;
      next_dn++;
    }
  }

  // sanity check
  if (next_up!=num_upholes_ && next_dn!=num_dnholes_) {
    throw std::logic_error(" BasisState::set_random: ivalid state status-3");
  }
  n_up=0;
  n_dn=0;
  for (const auto& state : site_states_) {
    if (state.have_uphole() && state.uphole_id()>=0) n_up++;
    if (state.have_dnhole() && state.dnhole_id()>=0) n_dn++;
  }
  if (n_up!=num_upholes_ && n_dn!=num_dnholes_) {
    throw std::logic_error(" BasisState::set_random: ivalid state status-4");
  }

  // in case there are no holes
  if (num_upholes_==0) {
    uphole_sites_.resize(1);
    uphole_sites_[0] = -1;
  }
  if (num_dnholes_==0) {
    dnhole_sites_.resize(1);
    dnhole_sites_[0] = -1;
  }
}

void BasisState::set_random_old(void)
{
  for (auto& state : site_states_) state.reset();
  std::vector<int> all_sites(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_sites[i] = i;
  std::shuffle(all_sites.begin(),all_sites.end(),rng_);
  // UP spins & holes
  for (int i=0; i<num_upspins_; ++i) {
    int site = all_sites[i];
    site_states_[site].put_upspin(i);
    upspin_sites_[i] = site;
  }
  int j = 0;
  for (int i=num_upspins_; i<num_sites_; ++i) {
    int site = all_sites[i];
    site_states_[site].put_uphole(j);
    uphole_sites_[j++] = site;
  }
  // DN spins & holes
  std::shuffle(all_sites.begin(),all_sites.end(),rng_);
  int id = 0;
  j = 0;
  for (int i=0; i<num_sites_; ++i) {
    int site = all_sites[i];
    if (id < num_dnspins_) {
      if (site_projection_[site]==pjn_t::DOUBLON && site_states_[site].have_upspin()) {
        site_states_[site].put_dnhole(j);
        dnhole_sites_[j++] = site;
      }
      else {
        site_states_[site].put_dnspin(id);
        dnspin_sites_[id++] = site;
      }
    }
    else {
      site_states_[site].put_dnhole(j);
      dnhole_sites_[j++] = site;
    }
  }
  // consistency check
  if (id != num_dnspins_) {
    throw std::logic_error("* BasisState::set_random: consistency check failed!");
  }
  if (j != num_dnholes_) {
    throw std::logic_error("* BasisState::set_random: consistency check failed!");
  }
  // number of doubly occupied sites
  /*num_dblocc_sites_ = 0;
  for (int i=0; i<num_sites_; ++i) {
    if (site_projection_[i] != proj_t::null) {
      if (site_states_[i].occupancy()==2) num_dblocc_sites_++;
    }
  }*/
  // in case there are no holes
  if (num_upholes_==0) {
    uphole_sites_.resize(1);
    upspin_sites_[0] = -1;
  }
  if (num_dnholes_==0) {
    dnhole_sites_.resize(1);
    dnspin_sites_[0] = -1;
  }
}


void BasisState::set_custom(void)
{
  for (auto& state : site_states_) state.reset();
  std::vector<int> all_sites(num_sites_);
  for (int i=0; i<num_sites_; ++i) all_sites[i] = i;
  int j=0;
  for (int i=0; i<num_sites_; i+=2) all_sites[j++] = i;
  for (int i=1; i<num_sites_; i+=2) all_sites[j++] = i;
  // UP spins & holes
  for (int i=0; i<num_upspins_; ++i) {
    int site = all_sites[i];
    site_states_[site].put_upspin(i);
    upspin_sites_[i] = site;
  }
  j = 0;
  for (int i=num_upspins_; i<num_sites_; ++i) {
    int site = all_sites[i];
    site_states_[site].put_uphole(j);
    uphole_sites_[j++] = site;
  }
  // DN spins & holes
  j = 0;
  for (int i=1; i<num_sites_; i+=2) all_sites[j++] = i;
  for (int i=0; i<num_sites_; i+=2) all_sites[j++] = i;
  int id = 0;
  j = 0;
  for (int i=0; i<num_sites_; ++i) {
    int site = all_sites[i];
    if (id<num_dnspins_ && site_projection_[site] != pjn_t::DOUBLON) {
      site_states_[site].put_dnspin(id);
      dnspin_sites_[id] = site;
      id++;
    }
    else {
      site_states_[site].put_dnhole(j);
      dnhole_sites_[j++] = site;
    }
  }
  // consistency check
  if (id != num_dnspins_) {
    throw std::logic_error("* BasisState::set_custom: consistency check failed!");
  }
  if (j != num_dnholes_) {
    throw std::logic_error("* BasisState::set_custom: consistency check failed!");
  }

  // in case there are no holes
  if (num_upholes_==0) {
    uphole_sites_.resize(1);
    uphole_sites_[0] = -1;
  }
  if (num_dnholes_==0) {
    dnhole_sites_.resize(1);
    dnhole_sites_[0] = -1;
  }
}

bool BasisState::gen_upspin_hop(void)
{
  if (num_upholes_==0 || num_upspins_==0) {
    proposed_move_ = move_t::null;
    return false;
  }

  // choose a spin
  mv_upspin_ = rng_.random_upspin();
  up_frsite_ = upspin_sites_[mv_upspin_]; 

  // choose a hole
  mv_uphole_ = rng_.random_uphole();
  up_tosite_ = uphole_sites_[mv_uphole_]; 

  // check projection constraint
  if (have_projection_) {
    if (site_projection_[up_frsite_]==pjn_t::HOLON) {
      if (site_states_[up_frsite_].occupancy()==1) {
        proposed_move_ = move_t::null;
        return false;
      }
    }
    if (site_projection_[up_tosite_]==pjn_t::DOUBLON) {
      if (site_states_[up_tosite_].occupancy()==1) {
        proposed_move_ = move_t::null;
        return false;
      }
    }
  }

  // move possible
  proposed_move_=move_t::upspin_hop;
  return true;
}

bool BasisState::gen_dnspin_hop(void)
{
  if (num_dnholes_==0 || num_dnspins_==0) {
    proposed_move_ = move_t::null;
    return false;
  }

  // choose a spin
  mv_dnspin_ = rng_.random_dnspin();
  dn_frsite_ = dnspin_sites_[mv_dnspin_]; 

  // choose a hole
  mv_dnhole_ = rng_.random_dnhole();
  dn_tosite_ = dnhole_sites_[mv_dnhole_]; 

  // check projection constraint
  if (have_projection_) {
    if (site_projection_[dn_frsite_]==pjn_t::HOLON) {
      if (site_states_[dn_frsite_].occupancy()==1) {
        proposed_move_ = move_t::null;
        return false;
      }
    }
    if (site_projection_[dn_tosite_]==pjn_t::DOUBLON) {
      if (site_states_[dn_tosite_].occupancy()==1) {
        proposed_move_ = move_t::null;
        return false;
      }
    }
  }

  // move possible
  proposed_move_=move_t::dnspin_hop;
  return true;
}

bool BasisState::gen_exchange_move(void)
{
  if (num_upholes_==0 || num_upspins_==0) return false;
  if (num_dnholes_==0 || num_dnspins_==0) return false;
  mv_upspin_ = rng_.random_upspin();
  mv_dnspin_ = rng_.random_dnspin();
  up_frsite_ = upspin_sites_[mv_upspin_]; 
  up_tosite_ = dnspin_sites_[mv_dnspin_];
  dn_frsite_ = dnspin_sites_[mv_dnspin_]; 
  dn_tosite_ = upspin_sites_[mv_upspin_];
  if (site_states_[up_tosite_].have_upspin()) {
    proposed_move_ = move_t::null;
    return false;
  }
  if (site_states_[dn_tosite_].have_dnspin()) {
    proposed_move_ = move_t::null;
    return false;
  }
  proposed_move_ = move_t::exchange;
  return true;
}

const int& BasisState::which_upspin(void) const
{
  if (proposed_move_==move_t::upspin_hop || proposed_move_==move_t::exchange) {
    return mv_upspin_;
  }
  else {
    throw std::logic_error("BasisState::which_upspin: no upspin move exists");
  }
}

const int& BasisState::which_dnspin(void) const
{
  if (proposed_move_==move_t::dnspin_hop || proposed_move_==move_t::exchange) {
    return mv_dnspin_;
  }
  else {
    throw std::logic_error("BasisState::which_dnspin: no dnspin move exists");
  }
}

const int& BasisState::which_frsite(void) const
{
  if (proposed_move_==move_t::upspin_hop) {
    return up_frsite_;
  }
  else if (proposed_move_==move_t::dnspin_hop) {
    return dn_frsite_;
  }
  else {
    throw std::logic_error("BasisState::which_frsite: no existing move");
  }
}

const int& BasisState::which_site(void) const
{
  if (proposed_move_==move_t::upspin_hop) {
    return up_tosite_;
  }
  else if (proposed_move_==move_t::dnspin_hop) {
    return dn_tosite_;
  }
  else {
    throw std::logic_error("BasisState::which_site: no existing move");
  }
}

const int& BasisState::which_upsite(void) const
{
  if (proposed_move_==move_t::upspin_hop || proposed_move_==move_t::exchange) {
    return up_tosite_;
  }
  else {
    throw std::logic_error("BasisState::which_upsite: no existing move");
  }
}

const int& BasisState::which_dnsite(void) const
{
  if (proposed_move_==move_t::dnspin_hop || proposed_move_==move_t::exchange) {
    return dn_tosite_;
  }
  else {
    throw std::logic_error("BasisState::which_dnsite: no existing move");
  }
}

void BasisState::accept_last_move(void)
{
  /*
  if (temporary_move_ != move_t::null) {
    throw std::range_error(" BasisState::accept_last_move: temporary change found, serious problem");
  }*/

  switch (proposed_move_) {
    case move_t::upspin_hop:
      //num_dblocc_sites_ += dblocc_increament_;
      site_states_[up_frsite_].put_uphole(mv_uphole_);
      uphole_sites_[mv_uphole_] = up_frsite_;
      site_states_[up_tosite_].put_upspin(mv_upspin_);
      upspin_sites_[mv_upspin_] = up_tosite_;
      proposed_move_ = move_t::null;
      break;
    case move_t::dnspin_hop:
      //num_dblocc_sites_ += dblocc_increament_;
      site_states_[dn_frsite_].put_dnhole(mv_dnhole_);
      dnhole_sites_[mv_dnhole_] = dn_frsite_;
      site_states_[dn_tosite_].put_dnspin(mv_dnspin_);
      dnspin_sites_[mv_dnspin_] = dn_tosite_;
      proposed_move_ = move_t::null;
      break;
    case move_t::exchange:
      mv_uphole_ = site_states_[up_tosite_].uphole_id();
      mv_dnhole_ = site_states_[dn_tosite_].dnhole_id();
      // spin moves
      site_states_[up_tosite_].put_upspin(mv_upspin_);
      upspin_sites_[mv_upspin_] = up_tosite_;
      site_states_[dn_tosite_].put_dnspin(mv_dnspin_);
      dnspin_sites_[mv_dnspin_] = dn_tosite_;
      // holes move
      site_states_[dn_tosite_].put_uphole(mv_uphole_);
      uphole_sites_[mv_uphole_] = dn_tosite_;
      site_states_[up_tosite_].put_dnhole(mv_dnhole_);
      dnhole_sites_[mv_dnhole_] = up_tosite_;
      proposed_move_ = move_t::null;
      break;
    case move_t::null:
      break;
  }
}

/*
void BasisState::temporary_upspin_hop(const int& fr_site, const int& to_site) const
{
  if (temporary_move_ != move_t::null) {
    throw std::range_error(" BasisState::temporary_upspin_hop: already contains temporary change");
  }
  tmp_up_frsite_ = fr_site;
  tmp_up_tosite_ = to_site;
  tmp_mv_upspin_ = site_states_[fr_site].upspin_id(); 
  tmp_mv_uphole_ = site_states_[to_site].uphole_id(); 
  site_states_[fr_site].put_uphole(tmp_mv_uphole_);
  site_states_[to_site].put_upspin(tmp_mv_upspin_);
  temporary_move_ = move_t::upspin_hop;
}

void BasisState::undo_upspin_hop(void) const
{
  if (temporary_move_ != move_t::upspin_hop) {
    throw std::range_error(" BasisState::undo_upspin_hop: no 'upspin_hop' to undo");
  }
  site_states_[tmp_up_frsite_].put_upspin(tmp_mv_upspin_);
  site_states_[tmp_up_tosite_].put_uphole(tmp_mv_uphole_);
  temporary_move_ = move_t::null;
}

void BasisState::temporary_dnspin_hop(const int& fr_site, const int& to_site) const
{
  if (temporary_move_ != move_t::null) {
    throw std::range_error(" BasisState::temporary_upspin_hop: already contains temporary change");
  }
  tmp_dn_frsite_ = fr_site;
  tmp_dn_tosite_ = to_site;
  tmp_mv_dnspin_ = site_states_[fr_site].dnspin_id(); 
  tmp_mv_dnhole_ = site_states_[to_site].dnhole_id(); 
  site_states_[fr_site].put_dnhole(tmp_mv_dnhole_);
  site_states_[to_site].put_dnspin(tmp_mv_dnspin_);
  temporary_move_ = move_t::dnspin_hop;
}

void BasisState::undo_dnspin_hop(void) const
{
  if (temporary_move_ != move_t::dnspin_hop) {
    throw std::range_error(" BasisState::undo_upspin_hop: no 'upspin_hop' to undo");
  }
  site_states_[tmp_dn_frsite_].put_dnspin(tmp_mv_dnspin_);
  site_states_[tmp_dn_tosite_].put_dnhole(tmp_mv_dnhole_);
  temporary_move_ = move_t::null;
}
*/

int BasisState::op_ni_up(const int& site) const
{
  return site_states_[site].num_upspins();
}

int BasisState::op_ni_dn(const int& site) const
{
  return site_states_[site].num_dnspins();
}

int BasisState::op_ni_updn(const int& site) const
{
  if (site_states_[site].occupancy() == 2) return 1;
  else return 0;
}

int BasisState::op_ni_dblon(const int& site) const
{
  if (site_states_[site].occupancy() == 2) return 1;
  else return 0;
}

int BasisState::op_ni_holon(const int& site) const
{
  if (site_states_[site].occupancy() == 0) return 1;
  else return 0;
}

int BasisState::op_Sz(const int& site) const
{
  return site_states_[site].num_upspins()-site_states_[site].num_dnspins();
}

/*-----------------------------------------------------------
* Apply: (c^\dag_{i\up} c_{j\up} + c^\dag_{j\up} c_{i\up})
*-----------------------------------------------------------*/
bool BasisState::op_cdagc_up_plus(const int& site_i, const int& site_j) const
{
  proposed_move_ = move_t::null;
  if (site_states_[site_i].have_uphole() && site_states_[site_j].have_upspin()) {
    up_frsite_ = site_j;
    up_tosite_ = site_i;
  }
  else if (site_states_[site_i].have_upspin() && site_states_[site_j].have_uphole()) {
    up_frsite_ = site_i;
    up_tosite_ = site_j;
  }
  else return false;

  // check projection constraint 
  if (have_projection_) {
    if (site_projection_[up_frsite_]==pjn_t::HOLON) {
      if (site_states_[up_frsite_].occupancy()==1) {
        return false;
      }
    }
    if (site_projection_[up_tosite_]==pjn_t::DOUBLON) {
      if (site_states_[up_tosite_].occupancy()==1) {
        return false;
      }
    }
  }

  // successfull hop
  mv_upspin_ = site_states_[up_frsite_].upspin_id(); 
  proposed_move_ = move_t::upspin_hop;
  return true;
}

/*-----------------------------------------------------------
* Apply: c^\dag_{i\up} c_{i\up}
*-----------------------------------------------------------*/
bool BasisState::op_cdagc_up(const int& fr_site, const int& to_site) const
{
  proposed_move_ = move_t::null;
  if (site_states_[fr_site].have_upspin() && site_states_[to_site].have_uphole()) {
    up_frsite_ = fr_site;
    up_tosite_ = to_site;
  }
  else return false;

  // check projection constraint 
  if (have_projection_) {
    if (site_projection_[up_frsite_]==pjn_t::HOLON) {
      if (site_states_[up_frsite_].occupancy()==1) {
        return false;
      }
    }
    if (site_projection_[up_tosite_]==pjn_t::DOUBLON) {
      if (site_states_[up_tosite_].occupancy()==1) {
        return false;
      }
    }
  }

  // successfull hop
  mv_upspin_ = site_states_[up_frsite_].upspin_id(); 
  proposed_move_ = move_t::upspin_hop;
  return true;
}


/*-----------------------------------------------------------
* Apply: (c^\dag_{i\dn} c_{j\dn} + c^\dag_{j\dn} c_{i\dn})
*-----------------------------------------------------------*/
bool BasisState::op_cdagc_dn_plus(const int& site_i, const int& site_j) const
{
  proposed_move_ = move_t::null;
  if (site_states_[site_i].have_dnhole() && site_states_[site_j].have_dnspin()) {
    dn_frsite_ = site_j;
    dn_tosite_ = site_i;
  }
  else if (site_states_[site_i].have_dnspin() && site_states_[site_j].have_dnhole()) {
    dn_frsite_ = site_i;
    dn_tosite_ = site_j;
  }
  else return false;

  // check projection constraint
  if (have_projection_) {
    if (site_projection_[dn_frsite_]==pjn_t::HOLON) {
      if (site_states_[dn_frsite_].occupancy()==1) {
        return false;
      }
    }
    if (site_projection_[dn_tosite_]==pjn_t::DOUBLON) {
      if (site_states_[dn_tosite_].occupancy()==1) {
        return false;
      }
    }
  }

  // successfull hop
  mv_dnspin_ = site_states_[dn_frsite_].dnspin_id(); 
  proposed_move_ = move_t::dnspin_hop;
  return true;
}


bool BasisState::op_cdagc_dn(const int& fr_site, const int& to_site) const
{
  proposed_move_ = move_t::null;
  if (site_states_[fr_site].have_dnspin() && site_states_[to_site].have_dnhole()) {
    dn_frsite_ = fr_site;
    dn_tosite_ = to_site;
  }
  else return false;

  // check projection constraint
  if (have_projection_) {
    if (site_projection_[dn_frsite_]==pjn_t::HOLON) {
      if (site_states_[dn_frsite_].occupancy()==1) {
        return false;
      }
    }
    if (site_projection_[dn_tosite_]==pjn_t::DOUBLON) {
      if (site_states_[dn_tosite_].occupancy()==1) {
        return false;
      }
    }
  }

  // successfull hop
  mv_dnspin_ = site_states_[dn_frsite_].dnspin_id(); 
  proposed_move_ = move_t::dnspin_hop;
  return true;
}

bool BasisState::op_exchange_ud(const int& site_i, const int& site_j) const
{
  proposed_move_ = move_t::null;
  const SiteState* state_i = &site_states_[site_i];
  const SiteState* state_j = &site_states_[site_j];
  // if any of the two sites doubly occupied, no exchange possible
  if (state_i->occupancy()==2 || state_j->occupancy()==2) return false;
  if (state_i->have_upspin() && state_j->have_dnspin()) {
    mv_upspin_ = state_i->upspin_id();
    up_tosite_ = site_j;
    mv_dnspin_ = state_j->dnspin_id();
    dn_tosite_ = site_i;
    proposed_move_ = move_t::exchange;
    return true;
  }
  else if (state_i->have_dnspin() && state_j->have_upspin()) {
    mv_upspin_ = state_j->upspin_id();
    up_tosite_ = site_i;
    mv_dnspin_ = state_i->dnspin_id();
    dn_tosite_ = site_j;
    proposed_move_ = move_t::exchange;
    return true;
  }
  else return false;
}

std::ostream& operator<<(std::ostream& os, const BasisState& bs)
{
  int len = 12;
  int i = 0;
  os << std::string(60,'-') << std::endl;
  os << "Basis state:\n";
  for (const auto& s : bs.site_states_) {
    os << "(" << s.state().to_string() << ") ";
    i++;
    if (i==len) {
      os << std::endl; i = 0;
    }
  }
  os << std::endl; 
  os << std::string(60,'-') << std::endl;
  return os;
}



} // end namespace vmc

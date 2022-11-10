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

BasisState::BasisState(void) 
{
  set_vaccuum(0);
} 

BasisState::BasisState(const int& num_sites) 
{
  set_vaccuum(num_sites);
} 

BasisState::BasisState(const int& num_sites, const bool& allow_dbl) 
{
  set_vaccuum(num_sites, allow_dbl);
} 

void BasisState::set_vaccuum(const int& num_sites, const bool& allow_dbl) 
{
  num_sites_ = num_sites;
  double_occupancy_ = allow_dbl;
  site_states_.resize(num_sites_);
  projections_.resize(num_sites_);
  site_projection_.resize(num_sites_);
  for (auto& p : site_projection_) p = pjn::NONE;


  if (double_occupancy_)  {
    for (auto& p : projections_) p = proj_t::partial;
  }
  else {
    for (auto& p : projections_) p = proj_t::full;
  } 
  upspin_sites_.clear();
  dnspin_sites_.clear();
  uphole_sites_.clear();
  dnhole_sites_.clear();
  // rng site generator
  if (num_sites_>0) rng_.set_site_generator(0,num_sites_-1);
}

void BasisState::clear(void)
{
  upspin_sites_.clear();
  dnspin_sites_.clear();
  uphole_sites_.clear();
  dnhole_sites_.clear();
} 

void BasisState::set_projection(const int& site, const proj_t& p) 
{ 
  projections_[site]=p; 
}

void BasisState::init_spins(const int& num_upspins, const int& num_dnspins)
{
  num_upspins_ = num_upspins;
  num_dnspins_ = num_dnspins;
  if (num_upspins_>num_sites_ || num_dnspins_>num_sites_)
    throw std::range_error("* BasisState::init_spins: spin number exceed capacity");
  int num_states = 0;
  for (int i=0; i<num_sites_; ++i) {
    if (projections_[i]==proj_t::full) num_states += 1;
    else num_states += 2;
  }
  if (num_upspins_+num_dnspins_>num_states) {
    throw std::range_error("* BasisState::init_spins: spin number exceed capacity");
  }
  num_upholes_ = num_sites_ - num_upspins;
  num_dnholes_ = num_sites_ - num_dnspins;
  // resizing
  upspin_sites_.resize(num_upspins_);
  dnspin_sites_.resize(num_dnspins_);
  uphole_sites_.resize(num_upholes_);
  dnhole_sites_.resize(num_dnholes_);
  // random generator
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

void BasisState::allow_double_occupancy(const bool& allow)
{
  double_occupancy_ = allow;
} 

void BasisState::set_random(void)
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
      if (projections_[site]==proj_t::full && site_states_[site].have_upspin()) {
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
  num_dblocc_sites_ = 0;
  for (int i=0; i<num_sites_; ++i) {
    if (projections_[i] != proj_t::null) {
      if (site_states_[i].occupancy()==2) num_dblocc_sites_++;
    }
  }
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
    if (id<num_dnspins_ && projections_[site] != proj_t::full) {
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
    throw std::logic_error("* BasisState::set_random: consistency check failed!");
  }
  if (j != num_dnholes_) {
    throw std::logic_error("* BasisState::set_random: consistency check failed!");
  }
  // number of doubly occupied sites
  num_dblocc_sites_ = 0;
  for (int i=0; i<num_sites_; ++i) {
    if (projections_[i] != proj_t::null) {
      if (site_states_[i].occupancy()==2) num_dblocc_sites_++;
    }
  }
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

bool BasisState::gen_upspin_hop(void)
{
  //if (proposed_move_!=move_t::null) undo_last_move();
  if (num_upholes_==0 || num_upspins_==0) {
    proposed_move_ = move_t::null;
    return false;
  }
  mv_upspin_ = rng_.random_upspin();
  mv_uphole_ = rng_.random_uphole();
  up_frsite_ = upspin_sites_[mv_upspin_]; 
  up_tosite_ = uphole_sites_[mv_uphole_]; 
  int n_occ = site_states_[up_tosite_].occupancy();
  if (projections_[up_tosite_]==proj_t::full && n_occ==1) {
    proposed_move_ = move_t::null;
    return false;
  }
  else {
    proposed_move_=move_t::upspin_hop;
    if (projections_[up_tosite_]!=proj_t::null) {
      dblocc_increament_ = n_occ; // must be 0 or 1
    }
    else dblocc_increament_ = 0;
    if (projections_[up_frsite_]!=proj_t::null) {
      if (site_states_[up_frsite_].have_dnspin()) dblocc_increament_--;
    }
    return true;
  }
}

bool BasisState::gen_dnspin_hop(void)
{
  //if (proposed_move_!=move_t::null) undo_last_move();
  if (num_dnholes_==0 || num_dnspins_==0) {
    proposed_move_ = move_t::null;
    return false;
  }
  mv_dnspin_ = rng_.random_dnspin();
  mv_dnhole_ = rng_.random_dnhole();
  dn_frsite_ = dnspin_sites_[mv_dnspin_]; 
  dn_tosite_ = dnhole_sites_[mv_dnhole_]; 
  int n_occ = site_states_[dn_tosite_].occupancy();
  if (projections_[dn_tosite_]==proj_t::full && n_occ==1) {
    proposed_move_ = move_t::null;
    return false;
  }
  else {
    proposed_move_=move_t::dnspin_hop;
    if (projections_[dn_tosite_]!=proj_t::null) {
      dblocc_increament_ = n_occ; // must be 0 or 1
    }
    else dblocc_increament_ = 0;
    if (projections_[dn_frsite_]!=proj_t::null) {
      if (site_states_[dn_frsite_].have_upspin()) dblocc_increament_--;
    }
    return true;
  }
}

bool BasisState::gen_exchange_move(void)
{
  //if (proposed_move_!=move_t::null) undo_last_move();
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
  dblocc_increament_ = 0;
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
    throw std::logic_error("BasisState::which_site: no existing move");
  }
}

const int& BasisState::which_dnsite(void) const
{
  if (proposed_move_==move_t::dnspin_hop || proposed_move_==move_t::exchange) {
    return dn_tosite_;
  }
  else {
    throw std::logic_error("BasisState::which_site: no existing move");
  }
}

void BasisState::accept_last_move(void)
{
  if (temporary_move_ != move_t::null) {
    throw std::range_error(" BasisState::accept_last_move: temporary change found, serious problem");
  }
  // double occupancy count
  switch (proposed_move_) {
    case move_t::upspin_hop:
      num_dblocc_sites_ += dblocc_increament_;
      site_states_[up_frsite_].put_uphole(mv_uphole_);
      uphole_sites_[mv_uphole_] = up_frsite_;
      site_states_[up_tosite_].put_upspin(mv_upspin_);
      upspin_sites_[mv_upspin_] = up_tosite_;
      proposed_move_ = move_t::null;
      break;
    case move_t::dnspin_hop:
      num_dblocc_sites_ += dblocc_increament_;
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

bool BasisState::op_cdagc_up_plus(const int& site_i, const int& site_j) const
{
  if (site_states_[site_i].have_uphole() && site_states_[site_j].have_upspin()) {
    up_frsite_ = site_j;
    up_tosite_ = site_i;
  }
  else if (site_states_[site_i].have_upspin() && site_states_[site_j].have_uphole()) {
    up_frsite_ = site_i;
    up_tosite_ = site_j;
  }
  else return false;
  int n_occ = site_states_[up_tosite_].occupancy();
  if (projections_[up_tosite_]==proj_t::full && n_occ==1) {
    return false;
  }
  mv_upspin_ = site_states_[up_frsite_].upspin_id(); 
  if (projections_[up_tosite_]!=proj_t::null) {
    dblocc_increament_ = n_occ; // must be 0 or 1
  }
  else dblocc_increament_ = 0;
  if (projections_[up_frsite_]!=proj_t::null) {
    if (site_states_[up_frsite_].have_dnspin()) dblocc_increament_--;
  }
  proposed_move_ = move_t::upspin_hop;
  return true;
}

bool BasisState::op_cdagc_up(const int& fr_site, const int& to_site) const
{
  if (site_states_[fr_site].have_upspin() && site_states_[to_site].have_uphole()) {
    up_frsite_ = fr_site;
    up_tosite_ = to_site;
  }
  else return false;
  int n_occ = site_states_[up_tosite_].occupancy();
  if (projections_[up_tosite_]==proj_t::full && n_occ==1) {
    return false;
  }
  mv_upspin_ = site_states_[up_frsite_].upspin_id(); 
  if (projections_[up_tosite_]!=proj_t::null) {
    dblocc_increament_ = n_occ; // must be 0 or 1
  }
  else dblocc_increament_ = 0;
  if (projections_[up_frsite_]!=proj_t::null) {
    if (site_states_[up_frsite_].have_dnspin()) dblocc_increament_--;
  }
  proposed_move_ = move_t::upspin_hop;
  return true;
}


/*
bool BasisState::op_cdagc_up(const int& frsite, const int& tosite) const
{
  if (!site_states_[frsite].have_upspin()) {
    return false;
  }
  if (frsite==tosite) {
    mv_upspin_ = site_states_[frsite].upspin_id(); 
    dblocc_increament_ = 0;
    proposed_move_ = move_t::upspin_hop;
    return true;
  }
  if (projections_[tosite]==proj_t::full && site_states_[tosite].have_dnspin()) {
    return false;
  }
  if (site_states_[tosite].have_uphole()) {
    mv_upspin_ = site_states_[frsite].upspin_id(); 
    int n_occ = site_states_[tosite].occupancy();
    if (projections_[tosite]!=proj_t::null) {
      dblocc_increament_ = n_occ; // must be 0 or 1
    }
    else dblocc_increament_ = 0;
    if (projections_[frsite]!=proj_t::null) {
      if (site_states_[frsite].have_dnspin()) dblocc_increament_--;
    }
    proposed_move_ = move_t::upspin_hop;
    return true;
  }
  return false;
} 
*/

bool BasisState::op_cdagc_dn_plus(const int& site_i, const int& site_j) const
{
  if (site_states_[site_i].have_dnhole() && site_states_[site_j].have_dnspin()) {
    dn_frsite_ = site_j;
    dn_tosite_ = site_i;
  }
  else if (site_states_[site_i].have_dnspin() && site_states_[site_j].have_dnhole()) {
    dn_frsite_ = site_i;
    dn_tosite_ = site_j;
  }
  else return false;
  int n_occ = site_states_[dn_tosite_].occupancy();
  if (projections_[dn_tosite_]==proj_t::full && n_occ==1) {
    return false;
  }
  mv_dnspin_ = site_states_[dn_frsite_].dnspin_id(); 
  if (projections_[dn_tosite_]!=proj_t::null) {
    dblocc_increament_ = n_occ; // must be 0 or 1
  }
  else dblocc_increament_ = 0;
  if (projections_[dn_frsite_]!=proj_t::null) {
    if (site_states_[dn_frsite_].have_upspin()) dblocc_increament_--;
  }
  proposed_move_ = move_t::dnspin_hop;
  return true;
}


bool BasisState::op_cdagc_dn(const int& fr_site, const int& to_site) const
{
  if (site_states_[fr_site].have_dnspin() && site_states_[to_site].have_dnhole()) {
    dn_frsite_ = fr_site;
    dn_tosite_ = to_site;
  }
  else return false;
  int n_occ = site_states_[dn_tosite_].occupancy();
  if (projections_[dn_tosite_]==proj_t::full && n_occ==1) {
    return false;
  }
  mv_dnspin_ = site_states_[dn_frsite_].dnspin_id(); 
  if (projections_[dn_tosite_]!=proj_t::null) {
    dblocc_increament_ = n_occ; // must be 0 or 1
  }
  else dblocc_increament_ = 0;
  if (projections_[dn_frsite_]!=proj_t::null) {
    if (site_states_[dn_frsite_].have_upspin()) dblocc_increament_--;
  }
  proposed_move_ = move_t::dnspin_hop;
  return true;
}


/*
bool BasisState::op_cdagc_dn(const int& frsite, const int& tosite) const
{
  if (!site_states_[frsite].have_dnspin()) {
    return false;
  }
  if (frsite==tosite) {
    mv_dnspin_ = site_states_[frsite].dnspin_id(); 
    dblocc_increament_ = 0;
    proposed_move_ = move_t::dnspin_hop;
    return true;
  }
  if (projections_[tosite]==proj_t::full && site_states_[tosite].have_upspin()) {
    return false;
  }
  if (site_states_[tosite].have_dnhole()) {
    mv_dnspin_ = site_states_[frsite].dnspin_id(); 
    int n_occ = site_states_[tosite].occupancy();
    if (projections_[tosite]!=proj_t::null) {
      dblocc_increament_ = n_occ; // must be 0 or 1
    }
    else dblocc_increament_ = 0;
    if (projections_[frsite]!=proj_t::null) {
      if (site_states_[frsite].have_upspin()) dblocc_increament_--;
    }
    proposed_move_ = move_t::dnspin_hop;
    return true;
  }
  return false;
}
*/

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

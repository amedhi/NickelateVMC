/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-02-13 10:16:02
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-12 21:26:24
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef BASISSTATE_H
#define BASISSTATE_H

#include <iostream>
#include <vector>
#include <utility>
#include <bitset>
#include <array>
#include "./random.h"
#include "../model/model.h"

namespace vmc {

enum class move_t {upspin_hop, dnspin_hop, exchange, null};

class SiteState 
{
public:
  SiteState() { reset(); }
  void reset(void) 
  {
    state_.reset();
    spin_id_.resize(2);
    spin_id_[UP] = -1;
    spin_id_[DN] = -1;
  }
  void put_upspin(const int& n) {state_.set(UP); spin_id_[UP]=n; }
  void put_dnspin(const int& n) {state_.set(DN); spin_id_[DN]=n; } 
  void put_uphole(const int& n) {state_.reset(UP); spin_id_[UP]=n;} 
  void put_dnhole(const int& n) {state_.reset(DN); spin_id_[DN]=n;} 
  const std::bitset<2>& state(void) const { return state_; }
  const int& upspin_id(void) const { return spin_id_[UP]; }
  const int& dnspin_id(void) const { return spin_id_[DN]; }
  const int& uphole_id(void) const { return spin_id_[UP]; }
  const int& dnhole_id(void) const { return spin_id_[DN]; }
  int num_upspins(void) const { return state_[UP]; }  
  int num_dnspins(void) const { return state_[DN]; }  
  bool have_upspin(void) const { return state_.test(UP); }  
  bool have_dnspin(void) const { return state_.test(DN); }  
  bool have_uphole(void) const { return !state_.test(UP); }  
  bool have_dnhole(void) const { return !state_.test(DN); }  
  bool is_full(void) const { return state_.all(); }  
  bool is_empty(void) const { return state_.none(); }  
  bool is_not_empty(void) const { return state_.any(); }  
  int occupancy(void) const { return state_.count(); }  
private:
  enum spin {UP, DN};
  std::bitset<2> state_;
  std::vector<int> spin_id_;
};

class BasisState 
{
public:
  BasisState() {} 
  BasisState(const lattice::Lattice& lattice, const model::Hamiltonian& model)
    { init(lattice, model); }
  ~BasisState() {} 
  void init(const lattice::Lattice& lattice, const model::Hamiltonian& model);
  //void allow_double_occupancy(const bool& allow);
  void init_spins(const int& num_upspins, const int& num_dnspins);
  int num_atoms(void) const { return atom_list_.size(); }
  const lattice::Lattice::Atom& atom(const int& atomid) const { return atom_list_[atomid]; }
  const int& atomid(const int& site) const { return atomid_list_[site]; }
  const lattice::Lattice::Atom& parent_atom(const int& site) const { int aid = atomid_list_[site]; return atom_list_[aid]; }
  void set_random(void);
  void set_random_old(void);
  void set_custom(void);
  //const bool& double_occupancy(void) const { return double_occupancy_; }
  bool gen_upspin_hop(void);
  bool gen_dnspin_hop(void);
  bool gen_exchange_move(void);
  const int& which_upspin(void) const;
  const int& which_dnspin(void) const;
  const int& which_frsite(void) const; 
  const int& which_site(void) const; 
  const int& which_upsite(void) const; 
  const int& which_dnsite(void) const; 
  //const int& dblocc_count(void) const { return num_dblocc_sites_; }
  //const int& dblocc_increament(void) const { return dblocc_increament_; }
  void accept_last_move(void);
  int op_ni(const int& site) const;
  int op_ni_up(const int& site) const;
  int op_ni_dn(const int& site) const;
  int op_ni_updn(const int& site) const;
  int op_ni_dblon(const int& site) const;
  int op_ni_holon(const int& site) const;
  int op_Sz(const int& site) const;
  bool op_cdagc_up_plus(const int& site_i, const int& site_j) const;
  bool op_cdagc_up(const int& frsite, const int& tosite) const;
  bool op_cdagc_dn_plus(const int& site_i, const int& site_j) const;
  bool op_cdagc_dn(const int& frsite, const int& tosite) const;
  bool op_exchange_ud(const int& fr_site, const int& to_site) const;
  const std::vector<int>& upspin_sites(void) const { return upspin_sites_; }
  const std::vector<int>& dnspin_sites(void) const { return dnspin_sites_; }
  void temporary_upspin_hop(const int& fr_site, const int& to_site) const;
  void temporary_dnspin_hop(const int& fr_site, const int& to_site) const;
  void undo_upspin_hop(void) const;
  void undo_dnspin_hop(void) const;
  RandomGenerator& rng(void) const { return rng_; }
  friend std::ostream& operator<<(std::ostream& os, const BasisState& bs);
private:
  enum class pjn_t {DOUBLON, HOLON, NONE};
  bool have_projection_{false};
  mutable RandomGenerator rng_;
  int num_sites_{0};
  int num_upspins_{0};
  int num_dnspins_{0};
  int num_upholes_{0};
  int num_dnholes_{0};
  //int num_dblocc_sites_{0};
  //bool double_occupancy_{true};
  std::vector<SiteState> site_states_;
  std::vector<pjn_t> site_projection_;
  std::vector<int> upspin_sites_;
  std::vector<int> dnspin_sites_;
  std::vector<int> uphole_sites_;
  std::vector<int> dnhole_sites_;
  /*------------------------------------------------------- 
   * Store the list of 'atom'-s. In single-orbital systems, 
   * an 'atom' is same as a 'site'. In multi-orbital systems,
   * an 'atom' store the id-s of sites (orbitals) that belongs
   * to it. The information is needed in case of inter-orbital
   * interaction.
   *-------------------------------------------------------*/
  std::vector<int> atomid_list_; // of sites (orbitals)
  std::vector<lattice::Lattice::Atom> atom_list_; 

  mutable move_t proposed_move_;
  //mutable int dblocc_increament_{0};
  mutable int mv_upspin_;
  mutable int up_frsite_;
  mutable int up_tosite_;
  mutable int mv_uphole_;
  mutable int mv_dnspin_;
  mutable int dn_frsite_;
  mutable int dn_tosite_;
  mutable int mv_dnhole_;

  // temporary change of state
  /*
  mutable move_t temporary_move_{move_t::null};
  mutable int tmp_mv_upspin_;
  mutable int tmp_up_frsite_;
  mutable int tmp_up_tosite_;
  mutable int tmp_mv_uphole_;
  mutable int tmp_mv_dnspin_;
  mutable int tmp_dn_frsite_;
  mutable int tmp_dn_tosite_;
  mutable int tmp_mv_dnhole_;
  */

  void clear(void); 
};




} // end namespace vmc

#endif

/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-19 23:06:41
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 11:32:25
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iomanip>
#include "./bcs_state.h"

namespace var {

BCS_State::BCS_State(const bcs& order_type, const input::Parameters& inputs, 
    const lattice::LatticeGraph& graph) 
  : GroundState(true)
{
  init(order_type, inputs, graph);
}

int BCS_State::init(const bcs& order_type, const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph)
{
  name_ = "BCS";
  // sites & bonds
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  // particle number
  set_particle_num(inputs);
  // infinity limit
  large_number_ = 5.0E+4;

  // build MF Hamiltonian
  order_type_ = order_type;
  varparms_.clear();
  mf_model_.init(graph.lattice());
  std::string name;
  double defval, lb, ub, dh;
  using namespace model;
  model::CouplingConstant cc;
  if (order_type_==bcs::swave) {
    order_name_ = "s-wave";

    if (graph.lattice().id()==lattice::lattice_id::SQUARE) {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
      // s-wave site pairing
      mf_model_.add_siteterm(name="pairing", cc="delta_sc", op::pair_create());
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      // variational parameters
      defval = mf_model_.get_parameter_value("delta_sc");
      varparms_.add("delta_sc",defval,lb=1.0E-4,ub=2.0,dh=0.01);
    }

    else if (graph.lattice().id()==lattice::lattice_id::NICKELATE_2L) {
      mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="tp", defval=0.0, inputs);
      mf_model_.add_parameter(name="th", defval=0.0, inputs);
      mf_model_.add_parameter(name="delta_N", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_R", defval=1.0, inputs);
      mf_model_.add_parameter(name="mu_N", defval=0.0, inputs);
      mf_model_.add_parameter(name="mu_R", defval=0.0, inputs);

      // bond operators
      cc.create(7);
      cc.add_type(0, "t");
      cc.add_type(1, "t");
      cc.add_type(2, "tp");
      cc.add_type(3, "t");
      cc.add_type(4, "t");
      cc.add_type(5, "tp");
      cc.add_type(6, "th");
      mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

      // site operators
      cc.create(2);
      cc.add_type(0, "-mu_N");
      cc.add_type(1, "e_R-mu_R");
      mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());

      // pairing term
      // s-wave paring on Ni & R sites
      cc.create(2);
      cc.add_type(0, "delta_N");
      cc.add_type(1, "delta_R");
      mf_model_.add_siteterm(name="singlet", cc, op::pair_create());

      // variational parameters
      defval = mf_model_.get_parameter_value("delta_N");
      varparms_.add("delta_N", defval,lb=1.0E-3,ub=3.0,dh=0.01);
      defval = mf_model_.get_parameter_value("delta_R");
      varparms_.add("delta_R",defval,lb=1.0E-3,ub=3.0,dh=0.01);
      defval = mf_model_.get_parameter_value("mu_N");
      varparms_.add("mu_N", defval,lb=-2.0,ub=+2.0,dh=0.1);
      defval = mf_model_.get_parameter_value("mu_R");
      varparms_.add("mu_R", defval,lb=-2.0,ub=+2.0,dh=0.1);
      noninteracting_mu_ = false;
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for this lattice");
    }
  }
  else if (order_type_==bcs::extended_swave) {
    order_name_ = "extended_s-wave";
    if (graph.lattice().id()==lattice::lattice_id::SW_GRAPHENE) {
      mf_model_.add_parameter(name="t0", defval=1.0, inputs);
      mf_model_.add_parameter(name="t1", defval=1.0, inputs);
      mf_model_.add_parameter(name="t2", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      // bond operator terms
      cc.create(3);
      for (int i=0; i<3; ++i) {
        cc.add_type(i, "-t"+std::to_string(i));
      }
      mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
      // pairing
      mf_model_.add_bondterm(name="pairing", cc="delta_sc", op::pair_create());
      // ch potemtial
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      // variational parameters
      defval = mf_model_.get_parameter_value("delta_sc");
      varparms_.add("delta_sc",defval,lb=1.0E-4,ub=2.0,dh=0.01);
    }

    else {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
      mf_model_.add_bondterm(name="pairing", cc="delta_sc", op::pair_create());
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      // variational parameters
      defval = mf_model_.get_parameter_value("delta_sc");
      varparms_.add("delta_sc",defval,lb=1.0E-4,ub=2.0,dh=0.01);
    }
  }
  else if (order_type_==bcs::dwave) {
    order_name_ = "d-wave";
    if (graph.lattice().id()==lattice::lattice_id::SQUARE_NNN) {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="tp", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
      mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      cc = CouplingConstant({0,"delta_sc"},{1,"-delta_sc"},{2,"0"},{3,"0"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
      double delta_sc = mf_model_.get_parameter_value("delta_sc");
      varparms_.add("delta_sc", defval=delta_sc, lb=0.0, ub=2.0, 0.01);
    }
    else if (graph.lattice().id()==lattice::lattice_id::SIMPLECUBIC) {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="tp", defval=1.0, inputs);
      mf_model_.add_parameter(name="tv", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"},{4,"-tv"});
      mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      cc = CouplingConstant({0,"delta_sc"},{1,"-delta_sc"},{2,"0"},{3,"0"},{4,"0"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
      double delta_sc = mf_model_.get_parameter_value("delta_sc");
      varparms_.add("delta_sc", defval=delta_sc, lb=0.0, ub=2.0,dh=0.01);
      double tv = mf_model_.get_parameter_value("tv");
      varparms_.add("tv", defval=tv, lb=-1.0, ub=1.0,dh=0.01);
    }
    else {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      cc = CouplingConstant({0, "delta_sc"}, {1, "-delta_sc"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
      // variational parameters
      varparms_.add("delta_sc", defval=1.0,lb=0.0,ub=2.0,dh=0.01);
    }
  }
  else if (order_type_==bcs::d_plus_id) {
    order_name_ = "d+id";
    if (graph.lattice().id()==lattice::lattice_id::SQUARE_NNN) {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      // pairing term
      double theta = 2.0*pi()/3.0;
      double theta2 = 2.0*theta;
      mf_model_.add_constant("theta",theta);
      mf_model_.add_constant("theta2",theta2);
      cc = CouplingConstant({0,"delta_sc"}, {1,"exp(i*theta)*delta_sc"}, {2,"exp(i*theta2)*delta_sc"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
      // variational parameters
      varparms_.add("delta_sc", defval=1.0, lb=0.0, ub=2.0);
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for this lattice");
    }
  }
  else if (order_type_==bcs::custom_sc) {
    order_name_ = "custom_sc";
    if (graph.lattice().id()==lattice::lattice_id::SQUARE_IONIC) {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="tp", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_d", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_xy", defval=0.0, inputs);
      mf_model_.add_parameter(name="P", defval=0.0, inputs);
      mf_model_.add_parameter(name="m", defval=0.0, inputs);

      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
      mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

      // SC pairing term
      cc = CouplingConstant({0,"delta_d"},{1,"-delta_d"},
        {2,"-delta_xy"},{3,"delta_xy"});
      //mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
      // site term
      cc.create(2);
      cc.add_type(0, "P-m");
      cc.add_type(1, "m-P");
      //mf_model_.add_siteterm(name="ni_up", cc, op::ni_up());
      cc.create(2);
      cc.add_type(0, "P+m");
      cc.add_type(1, "-(P+m)");
      //mf_model_.add_siteterm(name="ni_dn", cc, op::ni_dn());

      // variational parameters
      defval = mf_model_.get_parameter_value("delta_d");
      varparms_.add("delta_d",defval,lb=1.0E-3,ub=2.0,dh=0.01);
      defval = mf_model_.get_parameter_value("delta_xy");
      varparms_.add("delta_xy",defval,lb=1.0E-3,ub=2.0,dh=0.01);
      defval = mf_model_.get_parameter_value("P");
      varparms_.add("P", defval,lb=-1.0,ub=+1.0,dh=0.1);
      defval = mf_model_.get_parameter_value("m");
      varparms_.add("m", defval,lb=0.0,ub=+1.0,dh=0.1);
    }
    else if (graph.lattice().id()==lattice::lattice_id::NICKELATE) {
        mf_model_.add_parameter(name="e_N", defval=0.0, inputs);
        mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_NN_100", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_NN_001", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_NN_110", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_NN_200", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_NN_001", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_100", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_001", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_101", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_102", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_110", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_002", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_111", defval=0.0, inputs);
        mf_model_.add_parameter(name="tv_RN_200", defval=0.0, inputs);
        mf_model_.add_parameter(name="tv_RN_202", defval=0.0, inputs);
        mf_model_.add_parameter(name="delta_N", defval=1.0, inputs);
        mf_model_.add_parameter(name="delta_R", defval=1.0, inputs);
        mf_model_.add_parameter(name="mu_N", defval=0.0, inputs);
        mf_model_.add_parameter(name="mu_R", defval=0.0, inputs);

        // bond operators
        cc.create(15);
        cc.add_type(0, "t_NN_100");
        cc.add_type(1, "t_NN_100");
        cc.add_type(2, "t_NN_110");
        cc.add_type(3, "t_NN_200");
        cc.add_type(4, "t_NN_001");
        cc.add_type(5, "t_RR_100");
        cc.add_type(6, "t_RR_100");
        cc.add_type(7, "t_RR_001");
        cc.add_type(8, "t_RR_101");
        cc.add_type(9, "t_RR_102");
        cc.add_type(10, "t_RR_110");
        cc.add_type(11, "t_RR_002");
        cc.add_type(12, "t_RR_111");
        cc.add_type(13, "tv_RN_200");
        cc.add_type(14, "tv_RN_202");
        mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

        // site operators
        cc.create(2);
        cc.add_type(0, "e_N-mu_N");
        cc.add_type(1, "e_R-mu_R");
        mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());

        // pairing term
        // d-wave paring in Ni-layer
        cc.create(15);
        cc.add_type(0, "delta_N");
        cc.add_type(1, "-delta_N");
        cc.add_type(5, "delta_R");
        cc.add_type(6, "-delta_R");
        for (int i=2; i<5; ++i) cc.add_type(i, "0");
        for (int i=7; i<15; ++i) cc.add_type(i, "0");
        mf_model_.add_bondterm(name="bond_singlet", cc, op::pair_create());

        // local s-wave at R-sites
        /*
        cc.create(2);
        cc.add_type(0, "0");
        cc.add_type(1, "delta_s");
        mf_model_.add_siteterm(name="site_singlet", cc, op::pair_create());
        */

        // variational parameters
        defval = mf_model_.get_parameter_value("delta_N");
        varparms_.add("delta_N", defval,lb=1.0E-3,ub=2.0,dh=0.01);
        defval = mf_model_.get_parameter_value("delta_R");
        varparms_.add("delta_R",defval,lb=1.0E-3,ub=2.0,dh=0.01);
        defval = mf_model_.get_parameter_value("mu_N");
        varparms_.add("mu_N", defval,lb=-2.0,ub=+1.0,dh=0.1);
        defval = mf_model_.get_parameter_value("mu_R");
        varparms_.add("mu_R", defval,lb=-2.0,ub=+1.0,dh=0.1);
        defval = mf_model_.get_parameter_value("tv_RN_200");
        varparms_.add("tv_RN_200",defval,lb=-1.0,ub=+1.0,dh=0.01);
        defval = mf_model_.get_parameter_value("tv_RN_202");
        varparms_.add("tv_RN_202",defval,lb=-1.0,ub=+1.0,dh=0.01);
        noninteracting_mu_ = false;
    }
    else if (graph.lattice().id()==lattice::lattice_id::NICKELATE_2D) {
        mf_model_.add_parameter(name="e_N", defval=0.0, inputs);
        mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_NN_100", defval=1.0, inputs);
        mf_model_.add_parameter(name="t_NN_110", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RR_100", defval=1.0, inputs);
        mf_model_.add_parameter(name="t_RR_110", defval=0.0, inputs);
        mf_model_.add_parameter(name="t_RN_200", defval=0.0, inputs);
        mf_model_.add_parameter(name="delta_N", defval=1.0, inputs);
        mf_model_.add_parameter(name="delta_R", defval=1.0, inputs);
        mf_model_.add_parameter(name="mu_N", defval=0.0, inputs);
        mf_model_.add_parameter(name="mu_R", defval=0.0, inputs);

        // bond operators
	      cc.create(7);
        cc.add_type(0, "t_NN_100");
        cc.add_type(1, "t_NN_100");
        cc.add_type(2, "t_NN_110");
        cc.add_type(3, "t_RR_100");
        cc.add_type(4, "t_RR_100");
        cc.add_type(5, "t_RR_110");
        cc.add_type(6, "t_RN_200");
        mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

        // site operators
        cc.create(2);
        cc.add_type(0, "e_N-mu_N");
        cc.add_type(1, "e_R-mu_R");
        mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());

        // pairing term
        // d-wave paring in Ni-layer
        cc.create(7);
        cc.add_type(0, "delta_N");
        cc.add_type(1, "-delta_N");
        cc.add_type(2, "0");
        cc.add_type(3, "delta_R");
        cc.add_type(4, "-delta_R");
        cc.add_type(5, "0");
        cc.add_type(6, "0");
        mf_model_.add_bondterm(name="bond_singlet", cc, op::pair_create());

        // variational parameters
        defval = mf_model_.get_parameter_value("delta_N");
        varparms_.add("delta_N", defval,lb=1.0E-3,ub=2.0,dh=0.01);
        defval = mf_model_.get_parameter_value("delta_R");
        varparms_.add("delta_R",defval,lb=1.0E-3,ub=2.0,dh=0.01);
        defval = mf_model_.get_parameter_value("mu_N");
        varparms_.add("mu_N", defval,lb=-2.0,ub=+1.0,dh=0.1);
        defval = mf_model_.get_parameter_value("mu_R");
        varparms_.add("mu_R", defval,lb=-2.0,ub=+1.0,dh=0.1);
        noninteracting_mu_ = false;
    }
    else if (graph.lattice().id()==lattice::lattice_id::NICKELATE_2L) {
        mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
        mf_model_.add_parameter(name="t", defval=1.0, inputs);
        mf_model_.add_parameter(name="tp", defval=0.0, inputs);
        mf_model_.add_parameter(name="th", defval=0.0, inputs);
        mf_model_.add_parameter(name="delta_N", defval=1.0, inputs);
        mf_model_.add_parameter(name="delta_R", defval=1.0, inputs);
        mf_model_.add_parameter(name="mu_N", defval=0.0, inputs);
        mf_model_.add_parameter(name="mu_R", defval=0.0, inputs);

        // bond operators
        cc.create(7);
        cc.add_type(0, "t");
        cc.add_type(1, "t");
        cc.add_type(2, "tp");
        cc.add_type(3, "t");
        cc.add_type(4, "t");
        cc.add_type(5, "tp");
        cc.add_type(6, "th");
        mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

        // site operators
        cc.create(2);
        cc.add_type(0, "-mu_N");
        cc.add_type(1, "e_R-mu_R");
        mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());

        // pairing term
        // d-wave paring in Ni-layer
        cc.create(7);
        cc.add_type(0, "delta_N");
        cc.add_type(1, "-delta_N");
        cc.add_type(2, "0");
        cc.add_type(3, "delta_R");
        cc.add_type(4, "-delta_R");
        cc.add_type(5, "0");
        cc.add_type(6, "0");
        mf_model_.add_bondterm(name="bond_singlet", cc, op::pair_create());

        // variational parameters
        defval = mf_model_.get_parameter_value("delta_N");
        varparms_.add("delta_N", defval,lb=1.0E-3,ub=2.0,dh=0.01);
        defval = mf_model_.get_parameter_value("delta_R");
        varparms_.add("delta_R",defval,lb=1.0E-3,ub=2.0,dh=0.01);
        defval = mf_model_.get_parameter_value("mu_N");
        varparms_.add("mu_N", defval,lb=-2.0,ub=+2.0,dh=0.1);
        defval = mf_model_.get_parameter_value("mu_R");
        varparms_.add("mu_R", defval,lb=-2.0,ub=+2.0,dh=0.1);
        noninteracting_mu_ = false;
    }
    else if (graph.lattice().id()==lattice::lattice_id::HONEYCOMB) {
      mf_model_.add_parameter(name="t", defval=1.0, inputs);
      mf_model_.add_parameter(name="delta_sc", defval=1.0, inputs);
      mf_model_.add_parameter(name="theta1", defval=0.0, inputs);
      mf_model_.add_parameter(name="theta2", defval=0.0, inputs);
      // bond operators
      mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
      // site operators
      mf_model_.add_siteterm(name="mu_term", cc="-mu", op::ni_sigma());
      // pairing term
      cc = CouplingConstant({0,"delta_sc"}, {1,"delta_sc*exp(i*theta1)"}, 
        {2,"delta_sc*exp(i*theta2)"});
      mf_model_.add_bondterm(name="pairing", cc, op::pair_create());
      // variational parameters
      //defval = mf_model_.get_parameter_value("delta_sc");
      //varparms_.add("delta_sc", defval, lb=1.0E-3, ub=2.0, dh=0.01);
      defval = mf_model_.get_parameter_value("theta1");
      varparms_.add("theta1", defval, lb=0.0, ub=two_pi(), dh=0.2);
      defval = mf_model_.get_parameter_value("theta2");
      varparms_.add("theta2", defval, lb=0.0, ub=two_pi(), dh=0.2);
    }
    else {
      throw std::range_error("BCS_State::BCS_State: state undefined for this lattice");
    }
  }
  else {
    throw std::range_error("BCS_State::BCS_State: unidefined bcs order");
  }
  // chemical potential
  if (graph.lattice().id()==lattice::lattice_id::NICKELATE ||
      graph.lattice().id()==lattice::lattice_id::NICKELATE_2D ) {
    // do nothing
  }
  else {
    int info;
    mf_model_.add_parameter(name="mu", defval=0.0, inputs, info);
    if (info == 0) noninteracting_mu_ = false;
    else noninteracting_mu_ = true;
    if (inputs.set_value("mu_variational", false, info)) 
      varparms_.add("mu",defval=0.0,lb=-3.0,ub=+3.0,dh=0.1);
  }

  // finalize MF Hamiltonian
  mf_model_.finalize(graph);
  num_varparms_ = varparms_.size();

  // bloch basis
  blochbasis_.construct(graph);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(graph);

  //std::cout<< "mu_0 = " <<  get_noninteracting_mu() << "\n";

  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  delta_k_.resize(kblock_dim_,kblock_dim_);
  dphi_k_.resize(kblock_dim_,kblock_dim_);
  phi_k_.resize(num_kpoints_);
  work_k_.resize(num_kpoints_);
  for (unsigned k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
    work_k_[k].resize(kblock_dim_,kblock_dim_);
  } 
  return 0;
}

std::string BCS_State::info_str(void) const
{
  std::ostringstream info;
  info << "# Ground State: '"<<name_<<" ("<<order_name_<<")'\n";
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
    ++i;
  }
  mf_model_.update_terms();
}

void BCS_State::get_wf_amplitudes(Matrix& psi) 
{
  // k-space pair amplitudes
  if (kblock_dim_==1) {
    get_pair_amplitudes_oneband(phi_k_);
  }
  else {
    get_pair_amplitudes_multiband(phi_k_);
  }
  // 'lattice space' pair amplitudes
  get_pair_amplitudes_sitebasis(phi_k_, psi);
}

void BCS_State::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  int i=0; 
  for (const auto& p : varparms_) {
    double h = p.diff_h();
    double x = p.value();
    double x_fwd = x+h; 
    if (x_fwd > p.ubound()) x_fwd = p.ubound();
    mf_model_.update_parameter(p.name(), x_fwd);
    mf_model_.update_terms();
    if (kblock_dim_==1) get_pair_amplitudes_oneband(phi_k_);
    else get_pair_amplitudes_multiband(phi_k_);
    double x_bwd = x-h; 
    if (x_bwd < p.lbound()) x_bwd = p.lbound();
    mf_model_.update_parameter(p.name(), x_bwd);
    mf_model_.update_terms();
    if (kblock_dim_==1) get_pair_amplitudes_oneband(work_k_);
    else get_pair_amplitudes_multiband(work_k_);
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

void BCS_State::get_pair_amplitudes_multiband(std::vector<ComplexMatrix>& phi_k)
{
  // BCS pair amplitudes for multi-band system (INTRABAND pairing only)
  for (int k=0; k<num_kpoints_; ++k) {
    Vector3d kvec = blochbasis_.kvector(k);
    //-------------'+k block'-------------
    // hamiltonian in k-space
    mf_model_.construct_kspace_block(kvec);
    // diagonalize quadratic part
    es_k_up.compute(mf_model_.quadratic_spinup_block());

    // pairing part 
    delta_k_ = mf_model_.pairing_part();
    //-------------'-k block'-------------
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_up.eigenvectors().conjugate();
    //std::cout << "delta_k = " << delta_k_ << "\n"; getchar();
    // bcs ampitudes in rotated basis (assuming INTRABAND pairing only)
    dphi_k_.setZero();
    for (int i=0; i<kblock_dim_; ++i) {
      double ek = es_k_up.eigenvalues()[i] + es_minusk_up.eigenvalues()[i];
      double deltak_sq = std::norm(delta_k_(i,i));
      double ek_plus_Ek = ek + std::sqrt(ek*ek + 4.0*deltak_sq);
      /*
      if (ek<0.0) dphi_k_(i,i) = 1.0;
      else dphi_k_(i,i) = 0.0;
      */
      if (std::sqrt(deltak_sq)<1.0E-12 && ek<0.0) {
        dphi_k_(i,i) = large_number_ * std::exp(ii()*std::arg(delta_k_(i,i)));
        //std::cout << "large_number \n";
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
    phi_k[k] = work_ * es_minusk_up.eigenvectors().transpose();
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
    es_minusk_up.compute(mf_model_.quadratic_spinup_block());
    // assuming 'singlet pairing', see notes
    work_ = 0.5*(delta_k_ + mf_model_.pairing_part().transpose());
    //std::cout << work_ << "\n"; getchar();
    // transform pairing part
    delta_k_ = es_k_up.eigenvectors().adjoint() * work_ * 
      es_minusk_up.eigenvectors().conjugate();
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



} // end namespace var

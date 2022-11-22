/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Author: Amal Medhi
* Date:   2016-03-11 13:02:35
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-04-12 21:16:21
*----------------------------------------------------------------------------*/
#include <cmath>
#include "model.h"
#include <boost/algorithm/string.hpp>

namespace model {

int Hamiltonian::define_model(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  //int info;
  //unsigned ntypes;
  //std::vector<MatrixElement> matrix_elem(20);
  double defval;
  //unsigned sitetype, change, type, src_type, tgt_type;
  std::string name; //, matrixelem, op, qn, site, src, tgt, fact;
  //SiteBasis site_basis;
  //BasisDescriptor basis;
  //QuantumNumber::value_type min, max, step;
  CouplingConstant cc;

  // define the models 
  model_name = inputs.set_value("model", "HUBBARD");
  boost::to_upper(model_name);

  if (model_name == "HUBBARD") {
    mid = model_id::HUBBARD;
    if (lattice.id() == lattice::lattice_id::SQUARE_NNN) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }

    else if (lattice.id() == lattice::lattice_id::SQUARE_2SITE) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      //add_parameter(name="W", defval=0.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());

      // site operator terms
      //cc.create(2);
      //cc.add_type(0, "W");
      //cc.add_type(1, "-W");
      //add_siteterm(name="ni_sigma", cc, op::ni_sigma());
    }

    else if (lattice.id() == lattice::lattice_id::CHAIN_2SITE) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      cc = CouplingConstant({0,"-t"},{1,"-t"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }

    else if (lattice.id() == lattice::lattice_id::SQUARE_STRIPE) {
      add_parameter(name="tx", defval=1.0, inputs);
      add_parameter(name="ty1", defval=1.0, inputs);
      add_parameter(name="ty2", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);

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
      add_bondterm(name="hopping", cc, op::spin_hop());

      // Hubbard terms
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }

    else if (lattice.id() == lattice::lattice_id::SIMPLECUBIC) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      add_parameter(name="th", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"},{4,"-th"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }
    else if (lattice.id() == lattice::lattice_id::NICKELATE) {
      // model parameters
      add_parameter(name="U", defval=0.0, inputs);
      add_parameter(name="e_N", defval=0.0, inputs);
      add_parameter(name="e_R", defval=0.0, inputs);
      add_parameter(name="t_NN_100", defval=0.0, inputs);
      add_parameter(name="t_NN_001", defval=0.0, inputs);
      add_parameter(name="t_NN_110", defval=0.0, inputs);
      add_parameter(name="t_NN_200", defval=0.0, inputs);
      add_parameter(name="t_NN_001", defval=0.0, inputs);
      add_parameter(name="t_RR_100", defval=0.0, inputs);
      add_parameter(name="t_RR_001", defval=0.0, inputs);
      add_parameter(name="t_RR_101", defval=0.0, inputs);
      add_parameter(name="t_RR_102", defval=0.0, inputs);
      add_parameter(name="t_RR_110", defval=0.0, inputs);
      add_parameter(name="t_RR_002", defval=0.0, inputs);
      add_parameter(name="t_RR_111", defval=0.0, inputs);
      add_parameter(name="t_RN_200", defval=0.0, inputs);
      add_parameter(name="t_RN_202", defval=0.0, inputs);

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
      cc.add_type(13, "t_RN_200");
      cc.add_type(14, "t_RN_202");
      add_bondterm(name="hopping", cc, op::spin_hop());

      // site operators
      cc.create(2);
      cc.add_type(0, "e_N");
      cc.add_type(1, "e_R");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());

      // interaction
      cc.create(2);
      cc.add_type(0, "U");
      cc.add_type(1, "0");
      add_siteterm(name="hubbard", cc, op::hubbard_int());
    }
    else if (lattice.id() == lattice::lattice_id::NICKELATE_2D) {
      add_parameter(name="U", defval=0.0, inputs);
      add_parameter(name="e_N", defval=0.0, inputs);
      add_parameter(name="e_R", defval=0.0, inputs);
      add_parameter(name="t_NN_100", defval=1.0, inputs);
      add_parameter(name="t_NN_110", defval=0.0, inputs);
      add_parameter(name="t_RR_100", defval=1.0, inputs);
      add_parameter(name="t_RR_110", defval=0.0, inputs);
      add_parameter(name="t_RN_200", defval=0.0, inputs);
      // bond operators
      cc.create(7);
      cc.add_type(0, "t_NN_100");
      cc.add_type(1, "t_NN_100");
      cc.add_type(2, "t_NN_110");
      cc.add_type(3, "t_RR_100");
      cc.add_type(4, "t_RR_100");
      cc.add_type(5, "t_RR_110");
      cc.add_type(6, "t_RN_200");
      add_bondterm(name="hopping", cc, op::spin_hop());
      // site operators
      cc.create(2);
      cc.add_type(0, "e_N");
      cc.add_type(1, "e_R");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());
      // interaction
      cc.create(2);
      cc.add_type(0, "U");
      cc.add_type(1, "0");
      add_siteterm(name="hubbard", cc, op::hubbard_int());
    }

    else if (lattice.id() == lattice::lattice_id::NICKELATE_2L) {
      add_parameter(name="U", defval=0.0, inputs);
      add_parameter(name="e_R", defval=0.0, inputs);
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=0.0, inputs);
      add_parameter(name="th", defval=1.0, inputs);
      // bond operators
      cc.create(7);
      cc.add_type(0, "-t");
      cc.add_type(1, "-t");
      cc.add_type(2, "-tp");
      cc.add_type(3, "-t");
      cc.add_type(4, "-t");
      cc.add_type(5, "-tp");
      cc.add_type(6, "-th");
      add_bondterm(name="hopping", cc, op::spin_hop());
      // site operators
      cc.create(2);
      cc.add_type(0, "0");
      cc.add_type(1, "e_R");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());
      // interaction
      cc.create(2);
      cc.add_type(0, "U");
      cc.add_type(1, "0");
      add_siteterm(name="hubbard", cc, op::hubbard_int());
    }

    else if (lattice.id() == lattice::lattice_id::SW_GRAPHENE) {
      add_parameter(name="t0", defval=1.0, inputs);
      add_parameter(name="t1", defval=1.0, inputs);
      add_parameter(name="t2", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      cc.create(3);
      for (int i=0; i<3; ++i) {
        cc.add_type(i, "-t"+std::to_string(i));
      }
      add_bondterm(name="hopping", cc, op::spin_hop());
      // interaction
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }
    else {
      // model parameters
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      add_bondterm(name="hopping", cc="-t", op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }
  } 

  //------------------------HUBBARD IONIC-------------------------------------
  else if (model_name == "HUBBARD_IONIC") {
    mid = model_id::HUBBARD_IONIC;
    if (lattice.id() == lattice::lattice_id::SQUARE_4SITE) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      add_parameter(name="W", defval=0.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);
      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      add_siteterm(name="hubbard", cc="U", op::hubbard_int());

      // site operator terms
      cc.create(2);
      cc.add_type(0, "-0.5*W");
      cc.add_type(1, "0.5*W");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for this lattice");
    }
  }

  else if (model_name == "TJ_IONIC") {
    mid = model_id::TJ_IONIC;
    if (lattice.id() == lattice::lattice_id::SQUARE_4SITE) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      add_parameter(name="W", defval=0.0, inputs);
      add_parameter(name="U", defval=0.0, inputs);

      // projection operator
      ProjectionOp pjn;
      pjn.set({0,projection_t::HOLON}, {1,projection_t::DOUBLON});
      set_projection_op(pjn);

      // hopping term
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"});
      add_bondterm(name="hopping", cc, op::spin_hop());

      // exchange term
      cc.create(3);
      cc.add_type(0, "2.0*t*t/(U+W)");
      cc.add_type(1, "2.0*t*t/(U+W)");
      cc.add_type(2, "4.0*tp*tp/U");
      add_bondterm(name="exchange", cc, op::sisj_plus());

      // NN density-density terms
      cc.create(3);
      cc.add_type(0, "t*t/(U+W)-2.0*t*t/W");
      cc.add_type(1, "t*t/(U+W)-2.0*t*t/W");
      cc.add_type(2, "0");
      add_bondterm(name="ni_nj", cc, op::ni_nj());

      // Hubbard interaction
      cc.create(2);
      cc.add_type(0, "0.5*(U-W)");
      cc.add_type(1, "0.5*(U-W)");
      add_siteterm(name="hubbard", cc, op::hubbard_int());

      // extra onsite terms
      cc.create(2);
      cc.add_type(0, "8.0*tp*tp/U+2.0*t*t/W");
      cc.add_type(1, "6.0*t*t/W-0.5*(U-W)-2.0*t*t/(U+W)");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());

      // extra terms
      //add_siteterm(name="hubbard", cc="U", op::hubbard_int());
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for this lattice");
    }
  }

  //------------------------TJ-------------------------------------
  else if (model_name == "TJ") {
    mid = model_id::TJ;
    // projection operator
    //ProjectionOp pjn = projection_t::DOUBLON; 
    set_projection_op(projection_t::DOUBLON);
    if (lattice.id() == lattice::lattice_id::SQUARE) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="J", defval=0.0, inputs);
      // bond operator terms
      add_bondterm(name="hopping", cc="-t", op::spin_hop());
      add_bondterm(name="exchange", cc="J", op::sisj_plus());
    }
    else if (lattice.id() == lattice::lattice_id::SQUARE_NNN) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      add_parameter(name="J", defval=0.0, inputs);
      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      cc = CouplingConstant({0,"J"},{1,"J"},{2,"0"},{3,"0"});
      add_bondterm(name="exchange", cc="J", op::sisj_plus());
    }
    else if (lattice.id() == lattice::lattice_id::SIMPLECUBIC) {
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=1.0, inputs);
      add_parameter(name="th", defval=1.0, inputs);
      add_parameter(name="J", defval=0.0, inputs);
      // bond operator terms
      cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"},{4,"-th"});
      add_bondterm(name="hopping", cc, op::spin_hop());
      add_bondterm(name="exchange", cc="J", op::sisj_plus());
    }
    else if (lattice.id() == lattice::lattice_id::NICKELATE_2D) {
      add_parameter(name="J", defval=0.0, inputs);
      add_parameter(name="e_N", defval=0.0, inputs);
      add_parameter(name="e_R", defval=0.0, inputs);
      add_parameter(name="t_NN_100", defval=1.0, inputs);
      add_parameter(name="t_NN_110", defval=0.0, inputs);
      add_parameter(name="t_RR_100", defval=1.0, inputs);
      add_parameter(name="t_RR_110", defval=0.0, inputs);
      add_parameter(name="t_RN_200", defval=0.0, inputs);
      // bond operators
      cc.create(7);
      cc.add_type(0, "t_NN_100");
      cc.add_type(1, "t_NN_100");
      cc.add_type(2, "t_NN_110");
      cc.add_type(3, "t_RR_100");
      cc.add_type(4, "t_RR_100");
      cc.add_type(5, "t_RR_110");
      cc.add_type(6, "t_RN_200");
      add_bondterm(name="hopping", cc, op::spin_hop());
      // site operators
      cc.create(2);
      cc.add_type(0, "e_N");
      cc.add_type(1, "e_R");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());
      // interaction
      cc.create(7);
      cc.add_type(0, "J");
      cc.add_type(1, "J");
      cc.add_type(2, "0");
      cc.add_type(3, "0");
      cc.add_type(4, "0");
      cc.add_type(5, "0");
      cc.add_type(6, "0");
      add_bondterm(name="exchange", cc, op::sisj_plus());
    }
    else if (lattice.id() == lattice::lattice_id::NICKELATE_2L) {
      add_parameter(name="J", defval=0.0, inputs);
      add_parameter(name="e_R", defval=0.0, inputs);
      add_parameter(name="t", defval=1.0, inputs);
      add_parameter(name="tp", defval=0.0, inputs);
      add_parameter(name="th", defval=1.0, inputs);
      // bond operators
      cc.create(7);
      cc.add_type(0, "t");
      cc.add_type(1, "t");
      cc.add_type(2, "tp");
      cc.add_type(3, "t");
      cc.add_type(4, "t");
      cc.add_type(5, "tp");
      cc.add_type(6, "th");
      add_bondterm(name="hopping", cc, op::spin_hop());
      // site operators
      cc.create(2);
      cc.add_type(0, "0");
      cc.add_type(1, "e_R");
      add_siteterm(name="ni_sigma", cc, op::ni_sigma());
      // interaction
      cc.create(7);
      cc.add_type(0, "J");
      cc.add_type(1, "J");
      cc.add_type(2, "0");
      cc.add_type(3, "0");
      cc.add_type(4, "0");
      cc.add_type(5, "0");
      cc.add_type(6, "0");
      add_bondterm(name="exchange", cc, op::sisj_plus());
    }
    else {
      throw std::range_error("*error: modellibrary: model not defined for this lattice");
    }
  } // end model_name == "TJ"
  else {
    throw std::range_error("*error: modellibrary: undefined model");
  }

  // if the model has site disorder
  /*
  if (site_disorder) {
    add_disorder_term(name="disorder", op::ni_sigma());
  }*/
  
  return 0;
}

int Hamiltonian::construct(const input::Parameters& inputs, 
  const lattice::Lattice& lattice)
{
  init(lattice);
  define_model(inputs, lattice);
  finalize(lattice);
  return 0;
}


} // end namespace model

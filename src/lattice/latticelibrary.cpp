/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* All rights reserved.
* Date:   2016-01-17 21:32:15
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-06-19 21:37:15
*----------------------------------------------------------------------------*/
#include <stdexcept>
#include <iomanip>
#include <boost/algorithm/string.hpp>
#include "lattice.h"
//#include "graph.h"

namespace lattice {

// define the lattice
int Lattice::define_lattice(void) 
{
  using pos = Eigen::Vector3i;
  using vec = Eigen::Vector3d;
  unsigned type, ngb, src, tgt;
  vec a1, a2, a3, coord;
  pos src_offset, tgt_offset, cell;

  /*------------- 'SQUARE' lattice--------------*/
  if (lname == "SQUARE") {
    // type
    lid = lattice_id::SQUARE;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};

    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));

    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));

    // add bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
    add_bond(type=1, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
  }

  else if (lname == "SQUARE_NNN") {
    // type
    lid = lattice_id::SQUARE_NNN;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=3,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(-1,1,0));
  }

  else if (lname == "SQUARE_IONIC") {
    // type
    lid = lattice_id::SQUARE_IONIC;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    // basis vectors
    double a = 1.0;
    double b = std::sqrt(2.0)*a;
    set_basis_vectors(a1=vec(b,0,0), a2=vec(0,b,0), a3=vec(0,0,0));
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    add_basis_site(type=1, coord=vec(0.5*b,0.5*b,0));
    // add bonds
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
    add_bond(type=0,ngb=1,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=1,ngb=1,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,ngb=1,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));

    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=3,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=2,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=3,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
  }

  else if (lname == "SIMPLE_CUBIC") {
    // type
    lid = lattice_id::SIMPLECUBIC;
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,1));
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=3,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(-1,1,0));
    add_bond(type=4,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,0,1));
  }

  /*------------- 'CHAIN' lattice--------------*/
  else if (lname == "CHAIN") {
    lid = lattice_id::CHAIN;
    extent[dim2] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,0,0), a3=vec(0,0,0));
    // add sites
    add_basis_site(type=0, coord=vec(0,0,0));
    // add bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(1,0,0));
  }

  else if (lname == "HONEYCOMB") {
    // type
    lid = lattice_id::HONEYCOMB;
    // basis vectors
    double x = 0.5; 
    double y = -0.5*std::sqrt(3.0); 
    set_basis_vectors(a1=vec(x,y,0), a2=vec(x,-y,0), a3=vec(0,0,0));
    // sites
    add_basis_site(type=0, coord=vec(0,0,0));
    y = 1.0/std::sqrt(3.0); 
    add_basis_site(type=1, coord=vec(0,y,0));
    // bonds
    add_bond(type=0, ngb=1, src=0, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(0,0,0));
    add_bond(type=1, ngb=1, src=0, src_offset=pos(0,0,0), tgt=1, tgt_offset=pos(1,0,0));
    add_bond(type=2, ngb=1, src=1, src_offset=pos(0,0,0), tgt=0, tgt_offset=pos(0,1,0));
  }

  else if (lname == "NICKELATE") {
    lid = lattice_id::NICKELATE;
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,1));
    // sites
    add_basis_site(type=0, coord=vec(0,0,0));
    add_basis_site(type=1, coord=0.5*(basis_vector_a1()+
      basis_vector_a2()+basis_vector_a3()));

    // inter-cell Ni-Ni bonds
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));

    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,-1,0));

    add_bond(type=3,ngb=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(2,0,0));
    add_bond(type=3,ngb=3,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,2,0));

    add_bond(type=4,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,0,1));

    // inter-cell R-R bonds
    add_bond(type=5,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=6,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));

    add_bond(type=7,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,1));

    add_bond(type=8,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,1));
    add_bond(type=8,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,0,1));
    add_bond(type=8,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,1));
    add_bond(type=8,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,-1,1));

    add_bond(type=9,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,2));
    add_bond(type=9,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,0,2));
    add_bond(type=9,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,2));
    add_bond(type=9,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,-1,2));

    add_bond(type=10,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,1,0));
    add_bond(type=10,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,-1,0));

    add_bond(type=11,ngb=3,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,2));

    add_bond(type=12,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,1,1));
    add_bond(type=12,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,-1,1));
    add_bond(type=12,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,1,1));
    add_bond(type=12,ngb=4,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(-1,-1,1));

    // inter-cell Ni-R bonds
    add_bond(type=13,ngb=3,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(2,0,0));
    add_bond(type=13,ngb=3,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,2,0));

    add_bond(type=14,ngb=5,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(2,0,2));
    add_bond(type=14,ngb=5,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(-2,0,2));
    add_bond(type=14,ngb=5,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,2,2));
    add_bond(type=14,ngb=5,src=1,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,-2,2));
  }

  else if (lname == "NICKELATE_2D") {
    lid = lattice_id::NICKELATE_2D;
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));
    // sites
    add_basis_site(type=0, coord=vec(0,0,0));
    add_basis_site(type=1, coord=0.5*(basis_vector_a1()+basis_vector_a2()));

    // Ni-Ni bonds
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,-1,0));
    // R-R bonds
    add_bond(type=3,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=4,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
    add_bond(type=5,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,1,0));
    add_bond(type=5,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,-1,0));
    // Ni-R bonds
    add_bond(type=6,ngb=2,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=6,ngb=2,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
  }

  else if (lname == "NICKELATE_2L") {
    lid = lattice_id::NICKELATE_2L;
    // basis vectors
    set_basis_vectors(a1=vec(1,0,0), a2=vec(0,1,0), a3=vec(0,0,0));
    // sites
    add_basis_site(type=0, coord=vec(0,0,0));
    add_basis_site(type=1, coord=vec(0,0,1));

    // Ni-Ni bonds
    add_bond(type=0,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,0,0));
    add_bond(type=1,ngb=1,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(0,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,1,0));
    add_bond(type=2,ngb=2,src=0,src_offset=pos(0,0,0),tgt=0,tgt_offset=pos(1,-1,0));
    // R-R bonds
    add_bond(type=3,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,0,0));
    add_bond(type=4,ngb=1,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,1,0));
    add_bond(type=5,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,1,0));
    add_bond(type=5,ngb=2,src=1,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(1,-1,0));
    // Ni-R bonds
    add_bond(type=6,ngb=1,src=0,src_offset=pos(0,0,0),tgt=1,tgt_offset=pos(0,0,0));
  }

  else if (lname == "SW_GRAPHENE") {
    // type
    lid = lattice_id::SW_GRAPHENE;
    extent[dim3] = Extent{1, boundary_type::open, boundary_type::open, 0.0};
    // basis vectors
    double a = 1.0;
    double x = a*0.6034; 
    double y = a*0.7975; 
    set_basis_vectors(a1=vec(x,-y,0), a2=vec(x,y,0), a3=vec(0,0,0));
    // sites
    Vector3d a1 = basis_vector_a1();
    Vector3d a2 = basis_vector_a2();

    std::vector<int> C1(2);
    std::vector<int> C2(4);
    std::vector<int> C3(4);
    std::vector<int> C4(2);
    std::vector<int> C5(4);
    C1[0] = add_basis_site(type=0,coord=(0.564*a1+0.435*a2));
    C1[1] = add_basis_site(type=0,coord=(0.435*a1+0.564*a2));
    C2[0] = add_basis_site(type=0,coord=(0.497*a1+0.789*a2));
    C2[1] = add_basis_site(type=0,coord=(0.211*a1+0.497*a2));
    C2[2] = add_basis_site(type=0,coord=(0.789*a1+0.497*a2));
    C2[3] = add_basis_site(type=0,coord=(0.497*a1+0.211*a2));
    C3[0] = add_basis_site(type=0,coord=(0.701*a1+0.938*a2));
    C3[1] = add_basis_site(type=0,coord=(0.062*a1+0.299*a2));
    C3[2] = add_basis_site(type=0,coord=(0.938*a1+0.701*a2));
    C3[3] = add_basis_site(type=0,coord=(0.299*a1+0.062*a2));
    C4[0] = add_basis_site(type=0,coord=0.91*(a1+a2));
    C4[1] = add_basis_site(type=0,coord=0.09*(a1+a2));
    C5[0] = add_basis_site(type=0,coord=(0.318*a1+0.857*a2));
    C5[1] = add_basis_site(type=0,coord=(0.143*a1+0.682*a2));
    C5[2] = add_basis_site(type=0,coord=(0.857*a1+0.318*a2));
    C5[3] = add_basis_site(type=0,coord=(0.682*a1+0.143*a2));
    // bonds
    add_bond(type=0,src=C1[0],src_offset=pos(0,0,0),tgt=C1[1],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C1[1],src_offset=pos(0,0,0),tgt=C2[0],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C1[1],src_offset=pos(0,0,0),tgt=C2[1],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C1[0],src_offset=pos(0,0,0),tgt=C2[2],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C1[0],src_offset=pos(0,0,0),tgt=C2[3],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[0],src_offset=pos(0,0,0),tgt=C3[0],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[1],src_offset=pos(0,0,0),tgt=C3[1],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[2],src_offset=pos(0,0,0),tgt=C3[2],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[3],src_offset=pos(0,0,0),tgt=C3[3],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C3[0],src_offset=pos(0,0,0),tgt=C4[0],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C3[1],src_offset=pos(0,0,0),tgt=C4[1],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C3[2],src_offset=pos(0,0,0),tgt=C4[0],tgt_offset=pos(0,0,0));
    add_bond(type=1,src=C3[3],src_offset=pos(0,0,0),tgt=C4[1],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[0],src_offset=pos(0,0,0),tgt=C5[0],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[1],src_offset=pos(0,0,0),tgt=C5[1],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[2],src_offset=pos(0,0,0),tgt=C5[2],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C2[3],src_offset=pos(0,0,0),tgt=C5[3],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C5[0],src_offset=pos(0,0,0),tgt=C5[1],tgt_offset=pos(0,0,0));
    add_bond(type=2,src=C5[2],src_offset=pos(0,0,0),tgt=C5[3],tgt_offset=pos(0,0,0));
    // intercell
    add_bond(type=2,src=C3[0],src_offset=pos(0,0,0),tgt=C5[3],tgt_offset=pos(0,1,0));
    add_bond(type=2,src=C3[2],src_offset=pos(0,0,0),tgt=C5[1],tgt_offset=pos(1,0,0));
    add_bond(type=1,src=C4[0],src_offset=pos(0,0,0),tgt=C4[1],tgt_offset=pos(1,1,0));
    add_bond(type=2,src=C5[0],src_offset=pos(0,0,0),tgt=C3[3],tgt_offset=pos(0,1,0));
    add_bond(type=2,src=C5[2],src_offset=pos(0,0,0),tgt=C3[1],tgt_offset=pos(1,0,0));
  }


  /*------------- undefined lattice--------------*/
  else {
    throw std::range_error("error: latticelibrary: undefined lattice");
  }
  return 0;
}

// read lattice parameters
int Lattice::construct(const input::Parameters& parms) 
{

  int info;
  // name
  lname = parms.set_value("lattice", "NULL");
  boost::to_upper(lname);

  // sizes
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    std::string lsize = "lsize" + std::to_string(dim+1);
    extent[dim].size = parms.set_value(lsize, 1, info);
    if (extent[dim].size<1) throw std::range_error("error: latticelibrary: invalid lattice size");
  }

  // boundary conditions
  std::string bc; 
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    std::string lbc = "bc" + std::to_string(dim+1);
    bc = parms.set_value(lbc, "open", info);
    extent[dim].periodicity = boundary_condition(bc);
    extent[dim].bc = extent[dim].periodicity;
    if (extent[dim].bc == boundary_type::antiperiodic) extent[dim].bc = boundary_type::periodic;
    if (extent[dim].periodicity == boundary_type::antiperiodic) {
      extent[dim].bc_twist = pi();
    }
    else {
      extent[dim].bc_twist = 0.0;
    }
  }

  //extent[0].bc_twist = 0.0;
  //extent[1].bc_twist += two_pi()/extent[dim2].size * 2;
  //extent[2].bc_twist = 0.0;
  //for (unsigned dim=dim1; dim<=dim3; ++dim) {
  //  std::cout << "twist = " << extent[dim].bc_twist << "\n";
  //}
  //std::cout << "\n";

  
  // number of different boundary conditions
  int bc1_twists = parms.set_value("bc1_twists", 1, info);
  int bc2_twists = parms.set_value("bc2_twists", 1, info);
  int bc3_twists = parms.set_value("bc3_twists", 1, info);
  if (extent[dim1].size==1) bc1_twists = 1;
  if (extent[dim2].size==1) bc2_twists = 1;
  if (extent[dim3].size==1) bc3_twists = 1;
  num_total_twists_ = bc1_twists*bc2_twists*bc3_twists;
  twist_angles_.resize(num_total_twists_,3);
  // twist angles
  double dtheta1 = two_pi()/extent[dim1].size;
  double dtheta2 = two_pi()/extent[dim2].size;
  double dtheta3 = two_pi()/extent[dim3].size;
  int n = 0;
  for (int k=0; k<bc3_twists; ++k) {
    double step3 = k*dtheta3;
    for (int j=0; j<bc2_twists; ++j) {
      double step2 = j*dtheta2;
      for (int i=0; i<bc1_twists; ++i) {
        twist_angles_(n,0) = extent[dim1].bc_twist + i*dtheta1;
        twist_angles_(n,1) = extent[dim2].bc_twist + step2;
        twist_angles_(n,2) = extent[dim3].bc_twist + step3;
        n++;
      }
    }
  }
  // check
  //for (int i=0; i<num_total_twists_; ++i) {
  //  std::cout << "twist["<<i<<"] = "<<twist_angles_.row(i)<<"\n";
  //}
  //getchar();

  // empty unitcell
  Unitcell::clear();

  // impurities
  //impurity_sites_.clear();
  //impurity_bonds_.clear();

  define_lattice();
  finalize_lattice();
  //construct_graph();

  return 0;
}

int Lattice::finalize_lattice(void) 
{
  // Finalize the unit cell
  Unitcell::finalize();

  // copy the user set dimensions
  for (unsigned dim=dim1; dim<=dim3; ++dim) copy_extent[dim] = extent[dim];

  // is it necessary to construct 'symmetrized lattice'
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::open && extent[dim].size>1) {
      symmetrize_lattice();
      break;
    }
  }

  // number of unit cells & sites
  num_layer_cells_ = extent[dim1].size * extent[dim2].size;
  num_total_cells_ = num_layer_cells_ * extent[dim3].size;
  num_basis_sites_ = Unitcell::num_sites();
  num_total_sites_ = num_total_cells_ * num_basis_sites_;


  // check
  /*std::cout << "------Sites-------\n";
  for (unsigned i=0; i<unitcell.num_site(); ++i) {
    std::cout << Unitcell::site(i) << std::endl;
  }
  std::cout << "------Bonds-------\n";
  for (unsigned i=0; i<unitcell.num_bond(); ++i) {
    std::cout << Unitcell::bond(i) << std::endl;
  }*/

  // 'vector' & 'vector_id' attributes of the bonds
  std::map<int,unsigned> vecid_map;
  unsigned id=0;
  for (unsigned i=0; i<Unitcell::num_bonds(); ++i) {
    Vector3i ivec = Unitcell::bond(i).tgt().bravindex()-Unitcell::bond(i).src().bravindex();
    int key = ivec[0] + ivec[1]*extent[dim1].size + ivec[2]*num_layer_cells_;
    auto it = vecid_map.find(key);
    if (it != vecid_map.end()) Unitcell::bond(i).set_vector_id(it->second);
    else {
      vecid_map.insert({key, id});
      Unitcell::bond(i).set_vector_id(id);
      id++;
    }
    Unitcell::bond(i).set_vector((ivec[0]*vector_a1()+ivec[1]*vector_a2()+ivec[2]*vector_a3()));
    //std::cout << "bond " << i << ": vector_id = " << bond(i).vector_id() << "\n";
  }

  return 0;
}

int Lattice::symmetrize_lattice(void) 
{
  // initially, the 'dim' with periodic bc has size = 1
  spatial_dim = 0;
  Vector3d bvec;
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    if (extent[dim].bc==boundary_type::periodic) {
      spatial_dim++;
      // temporarily set size = 1 for dim with PBC
      extent[dim].size = 1;
      switch (dim) {
        case dim1: bvec = basis_vector_a1(); break;
        case dim2: bvec = basis_vector_a2(); break;
        case dim3: bvec = basis_vector_a3(); break;
      }
    }
  }
  // if 1 dimensional lattice, rotate the lattice to align 'bvec' along x-direction
  if (spatial_dim == 1) {
    // rotation matrix to do that
    Eigen::Matrix3d matrix = rotation_matrix(bvec, Vector3d(1.0, 0.0, 0.0));
    // rotate the unitcell
    Unitcell::rotate_by(matrix);
  }

  // number of unit cells & sites
  num_layer_cells_ = extent[dim1].size * extent[dim2].size;
  num_total_cells_ = num_layer_cells_ * extent[dim3].size;
  num_basis_sites_ = Unitcell::num_sites();
  num_total_sites_ = num_total_cells_ * num_basis_sites_;

  // Add the sites & the bonds to the symmetrized unitcell
  std::vector<Site> sites;
  std::vector<Bond> bonds;
  Unitcell translated_cell;
  Vector3i bravindex(0,0,0);
  for (unsigned i=0; i<num_total_cells_; ++i) {
    translated_cell = get_translated_cell(bravindex);
    // collect the sites
    for (unsigned n=0; n<translated_cell.num_sites(); ++n) 
      sites.push_back(translated_cell.site(n));
    // collect the bonds
    for (unsigned n=0; n<translated_cell.num_bonds(); ++n) 
      bonds.push_back(translated_cell.bond(n));
    bravindex = get_next_bravindex(bravindex);
  }

  // replace the old sites & bonds
  Unitcell::clear_sites();
  unsigned i = 0;
  for (auto& s : sites) {
    s.reset_uid(i++); 
    s.reset_bravindex(Vector3i(0,0,0));
    s.reset_cell_coord(Vector3d(0.0,0.0,0.0));
    Unitcell::add_site(s);
  }
  Unitcell::clear_bonds();
  for (auto& b : bonds) {
    if (connect_bond2(b, sites)) {
      b.reset_bravindex(Vector3i(0,0,0));
      Unitcell::add_bond(b);
    }
  }


  //std::cout << unitcell.vector_a1() << "\n";
  //std::cout << unitcell.vector_a2() << "\n";
  //std::cout << unitcell.vector_a3() << "\n";

  // extent & basis vectors of the symmetrized lattice
  for (unsigned dim=dim1; dim<=dim3; ++dim) {
    extent[dim] = copy_extent[dim];
    if (extent[dim].bc != boundary_type::periodic) {
      extent[dim].size = 1;
      switch (dim) {
        case dim1: Unitcell::reset_a1(Vector3d(0.0,0.0,0.0)); break;
        case dim2: Unitcell::reset_a2(Vector3d(0.0,0.0,0.0)); break;
        case dim3: Unitcell::reset_a3(Vector3d(0.0,0.0,0.0)); break;
      }
    }
  }
  //std::cout << unitcell.vector_a1() << "\n";
  //std::cout << unitcell.vector_a2() << "\n";
  //std::cout << unitcell.vector_a3() << "\n";

  return 0;
}






} // end namespace lattice

/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-19 16:03:53
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-19 23:14:39
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include <iomanip>
#include <string>
#include <stdexcept>
#include <Eigen/SVD>
#include <Eigen/QR>
#include <boost/algorithm/string.hpp>
#include <boost/filesystem.hpp>
#include "./optimizer.h"
//#include "./LBFGS/LBFGSB.h"
//#include "./LBFGS/LBFGSpp/Param.h"

namespace opt {

int Optimizer::init(const input::Parameters& inputs, const vmc::VMC& vmc) 
{
  // problem size
  num_parms_ = vmc.num_varp();
  vparms_.resize(num_parms_);
  vparms_start_.resize(num_parms_);
  lbound_ = vmc.varp_lbound();
  ubound_ = vmc.varp_ubound();
  range_ = ubound_-lbound_;
  grad_.resize(num_parms_);
  sr_matrix_.resize(num_parms_,num_parms_);

  return 0;
}

int Optimizer::optimize(vmc::VMC& vmc)
{
  //LBFGSpp::LBFGSBParam<double> param;
  //LBFGSpp::LBFGSBSolver<double> solver(param);


	return 0;
}




} // end name space optimizer

/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-19 07:07:40
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-20 11:12:26
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef OPTIMIZER_H
#define OPTIMIZER_H

#include <iostream>
#include <fstream>
#include <iomanip>
#include <deque>
//#include "./LBFGS/LBFGSB.h"
//#include "./LBFGS/LBFGSpp/Param.h"
#include "../vmc/vmc.h"


namespace opt {

class Optimizer
{
public:
  Optimizer() {} 
  Optimizer(const input::Parameters& parms); 
  ~Optimizer() {}
  int init(const input::Parameters& parms, const vmc::VMC& vmc);
  int optimize(vmc::VMC& vmc);
  const var::parm_vector& optimal_parms(void) const { return vparms_; }
  //const var::parm_vector& vp(void) { return varparms; }
private:
  int num_parms_;
  mcdata::MC_Observable optimal_parms_;
  mcdata::MC_Observable optimal_energy_;
  mcdata::MC_Observable energy_error_bar_;
  var::parm_vector vparms_;
  var::parm_vector vparms_start_;
  var::parm_vector lbound_;
  var::parm_vector ubound_;
  var::parm_vector range_;
  Eigen::MatrixXd sr_matrix_;
  Eigen::VectorXd grad_;
  std::vector<double> xvar_values_;
};


} // end namespace opt

#endif
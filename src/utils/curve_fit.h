/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2020-04-13 11:48:05
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2020-04-13 17:48:17
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef CURVEFIT_H
#define CURVEFIT_H

#include <iostream>
#include <cassert>
#include <Eigen/Dense>

namespace util {

using VectorXd = Eigen::VectorXd;
using MatrixXd = Eigen::MatrixXd;

/*
   "The Levenberg-Marquardt Method" of nonlinear curve fitting.
   See, Henri P. Gavin, "The Levenberg-Marquardt algorithm for 
   non-linear least square curve-fitting problems".
*/
class CurveFit
{
public:
	CurveFit() {}
	~CurveFit() {}
  template <typename Foo>
  inline int lmfit(Foo& func, VectorXd& p, const VectorXd& xdata, const VectorXd& ydata) 
  {
    assert(xdata.size()==ydata.size());
    VectorXd yerr = VectorXd::Ones(xdata.size());
    VectorXd perr = VectorXd::Zero(p.size());
    have_yerr = false;
    return lmfit(func, p, xdata, ydata, yerr, perr);
  }

	template <typename Foo>
  inline int lmfit(Foo& func, VectorXd& p, const VectorXd& xdata, const VectorXd& ydata, 
  	const VectorXd& yerr)
  {
  	assert(xdata.size()==ydata.size());
    assert(xdata.size()==yerr.size());
    VectorXd perr = VectorXd::Zero(p.size());
  	have_yerr = true;
  	return lmfit(func, p, xdata, ydata, yerr, perr);
  }

	template <typename Foo>
  inline int lmfit(Foo& func, VectorXd& p, const VectorXd& xdata, const VectorXd& ydata, 
  	const VectorXd& yerr, VectorXd& perr)
  {
  	assert(xdata.size()==ydata.size());
  	assert(xdata.size()==yerr.size());
  	assert(p.size()==perr.size());

  	m_data = xdata.size();
  	n_parm = p.size();
  	VectorXd weights = VectorXd::Ones(m_data);

  	MatrixXd JT(n_parm,m_data);
  	MatrixXd A(n_parm,m_data);
  	MatrixXd B(n_parm,n_parm);
  	VectorXd yfit(m_data);
  	VectorXd Y(n_parm);
  	VectorXd h_step(n_parm);
  	VectorXd p_next(n_parm);
    VectorXd hv(n_parm);
    VectorXd BD(n_parm);

  	double chisq = chi_square(func,p,xdata,ydata,weights);
    //std::cout << "chisq =" << chisq << "\n";
  	analytic_JacobianT(func,p,xdata,JT);
  	// Matrix A = J^T W 
    for (int i=0; i<m_data; ++i) A.col(i) = JT.col(i)*weights(i);
  	// Matrix J^T W J
   	B = A*JT.transpose();
    BD = B.diagonal();
  	// rhs vector in LM step
  	for (int i=0; i<m_data; ++i) yfit(i) = func(xdata[i],p);
    Y = A*(ydata-yfit);
    // iterations
    bool converged = false;
  	for (int iter=0; iter<200; ++iter) {
      // normalize B
  		for (int i=0; i<n_parm; ++i) B(i,i) *=  (1.0+lambda);

      h_step = B.fullPivLu().solve(Y);
      p_next = p + h_step;
      double chisq_next = chi_square(func,p_next,xdata,ydata,weights);
    	//std::cout << "chisq =" << chisq_next << "\n";
    	// LM metric
    	for (int i=0; i<n_parm; ++i) hv(i) = lambda*BD(i)*h_step(i)+Y(i);
    	double d = h_step.transpose() * hv;
    	double lm_rho = (chisq-chisq_next)/d;
    	//std::cout << "rho =" << lm_rho << "\n";
    	if (lm_rho > eps4) {
  			lambda = std::max(lambda/lambda_DN,10E-07);
    		p = p_next;
    		chisq = chisq_next;
  			analytic_JacobianT(func,p,xdata,JT);
    		for (int i=0; i<m_data; ++i) A.col(i) = JT.col(i)*weights(i);
   			B = A*JT.transpose();
    		BD = B.diagonal();
  			// rhs vector in LM step
  			for (int i=0; i<m_data; ++i) yfit(i) = func(xdata[i],p);
    		Y = A*(ydata-yfit);
    	}
    	else {
    		lambda = std::min(lambda*lambda_UP,1.0E+07);
    	}
    	// convergence
    	if (Y.lpNorm<Eigen::Infinity>() < eps1) {
    		converged = true;
    		break;
    	} 
    	double g = 0.0;
    	for (int i=0; i<n_parm; ++i) {
    		double r = std::abs(h_step(i)/p(i));
    		if (r > g) g = r;
    	}
    	if (g < eps2) {
    		converged = true;
    		break;
    	}
    	if (chisq/(m_data-n_parm+1) < eps3) {
    		converged = true;
    		break;
    	}
  	}
  	//std::cout << "p = " << p.transpose() << "\n";
  	// error estimate
   	B = A*JT.transpose();
   	MatrixXd covmat = B.inverse();
  	for (int i=0; i<n_parm; ++i) {
  		perr(i) = std::sqrt(covmat(i,i));
  	}
  	//std::cout << "perr = " << perr.transpose() << "\n";

  	return converged;
  }

private:
	int m_data;
	int n_parm;
	bool have_yerr{true};
  double lambda{1.0E-2};
  double lambda_UP{11};
  double lambda_DN{9};
  double eps1{1.0E-4};
  double eps2{1.0E-4};
  double eps3{1.0E-4};
  double eps4{1.0E-1};

	template <typename Foo>
  inline double chi_square(Foo& func, VectorXd& p, const VectorXd& xdata, 
  	const VectorXd& ydata, const VectorXd& weights) 
  {
  	double sum = 0.0;
  	for (int i=0; i<xdata.size(); ++i) {
  		double err = weights[i]*(ydata[i]-func(xdata[i],p));
  		sum += err*err;
  	}
  	return sum;
  }

	template <typename Foo>
  inline void analytic_JacobianT(Foo& func, VectorXd& p, const VectorXd& xdata, 
  	MatrixXd& JT)
  {
  	VectorXd pd(n_parm);
  	for (int i=0; i<m_data; ++i) {
  		func.derivative(xdata[i],p,pd);
  		JT.col(i) = pd;
  	}
  }
	
};

}  // namespace util
#endif
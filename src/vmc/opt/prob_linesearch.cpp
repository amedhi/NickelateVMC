/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-03-09 15:19:43
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-21 17:40:34
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <iostream>
#include <limits>
#include <string>
#include <random>
#include <stdexcept>
#include <Eigen/QR>
#include <boost/algorithm/string.hpp>
#include "./prob_linesearch.h"

namespace opt {
/*
Implementation of Probabilistic Line Search for Stochastic Optimization [1].
[1] M. Mahsereci and P. Hennig. Probabilistic line searches for stochastic
optimization. In Advances in Neural Information Processing Systems 28, pages
181-189, 2015.
*/

/*
The implementation is ported from the Python implementation of the same 
in https://github.com/ProbabilisticNumerics/probabilistic_line_search

Given below is the notice as required by the Liense for the software.
-----------------------------------------------------------------------
Copyright 2015, Max Planck Society.

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
-----------------------------------------------------------------------
*/


int ProbLineSearch::test(void) 
{
  std::cout << "Test CDF\n";
  std::cout << "CDF(0.0)=0.5 ? "<<CDF(0.0)<< "\n";
  std::cout << "CDF(1.0)=0.8413 ? "<<CDF(1.0) << "\n";
  std::cout << "CDF(3.0)=0.9987 ? "<<CDF(3.0) << "\n";
  std::cout << "CDF(-1.0)=0.1587 ? "<<CDF(-1.0) << "\n";
  std::cout << "CDF(-0.1)=0.4602 ? "<<CDF(-0.1) << "\n";
  std::cout << "\n";

  std::cout << "Test 'unbounded_bivariate_normal_integral'\n";
  double x;
  x = unbounded_bivariate_normal_integral(0.,0.,0.);
  std::cout << "ubni(0,0,0)=0.25 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(1.,0.,0.);
  std::cout << "ubni(1,0,0)=0.5 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(0.43, 2.5, -1.0);
  std::cout << "ubni(0.43, 2.5, -1.0)=0.0062 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(-0.17, 0.5, 1.0);
  std::cout << "ubni(-0.17, 0.5, 1.0)=0.0351 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(0., 0.5, 1.0);
  std::cout << "ubni(0., 0.5, 1.0)=0.0490 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(-1.0, 0., -3.0);
  std::cout << "ubni(-1.0, 0., -3.)=0.4987 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(0., -5., -5.);
  std::cout << "ubni(0., -5., -5.)=1.0 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(0., 5., 5.);
  std::cout << "ubni(0., 5., 5.)=0.0 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(1., 3., 3.);
  std::cout << "ubni(1., 3., 3.)=0.0013 ? "<<x<< "\n";
  x = unbounded_bivariate_normal_integral(-1., 0., 0.);
  std::cout << "ubni(-1., 0., 0.)=0.0 ? "<<x<< "\n";
  std::cout << "\n";

  std::cout << "Test 'bounded_bivariate_normal_integral'\n";
  x = bounded_bivariate_normal_integral(0.25, 0., 2.5, -1.2, 0.1);
  std::cout << "bvn(0.25, 0., 2.5, -1.2, 0.1)=0.1901 ? "<<x<< "\n";
  x = bounded_bivariate_normal_integral(0., 0., 1., 0., 1.);
  std::cout << "bvn(0., 0., 1., 0., 1.)=0.1165 ? "<<x<< "\n";
  x = bounded_bivariate_normal_integral(0.5, 0., 1., 0., 1.);
  std::cout << "bvn(0.5, 0., 1., 0., 1.)=0.1411 ? "<<x<< "\n";
  double inf = std::numeric_limits<double>::infinity();
  x = bounded_bivariate_normal_integral(0.5, 0., inf, 0., 1.);
  std::cout << "bvn(0.5, 0., inf, 0., 1.)=0.2059 ? "<<x<< "\n";

  return 0;
}

ProbLineSearch::ProbLineSearch(void)
{
  c1_ = 0.05;
  cW_ = 0.3; 
  fpush_ = 1.0; 
  alpha0_ = 0.02; // 0.01; 
  alpha_stats_ = alpha0_;
  target_df_ = 0.5;
  df_lo_ = -0.1, 
  df_hi_ = 1.1;
  max_steps_ = 10; 
  max_expl_ = 6;
  max_dmu0_ = 0.0; 
  max_change_factor_ = 10.0;
  expl_policy_ = "linear";
  //expl_policy_ = "exponential";

  num_steps_ = 0;
  num_expl_ = 0;
  new_search_ = false;
}

void ProbLineSearch::set_parameters(const double& c1, const double& cW, 
  const double& fpush, const double& alpha0, const double& target_df, 
  const double& df_lo, const double& df_hi, const int& max_steps, 
  const int& max_expl, const double& max_dmu0, 
  const double& max_change_factor, const std::string& expl_policy) 
{
  c1_ = c1; cW_=cW; fpush_=fpush; alpha0_=alpha0; target_df_=target_df;
  df_lo_=df_lo; df_hi_=df_hi; max_steps_=max_steps; max_expl_=max_expl;
  max_dmu0_=max_dmu0; max_change_factor_=max_change_factor;
  expl_policy_=expl_policy;
} 

double ProbLineSearch::start(const double& en, const double& en_err, 
  const RealVector& grad, const RealVector& grad_err, 
  const RealVector& search_dir)
{
  GP_.reset();
  num_steps_ = 0;
  num_expl_ = 0;
  abort_status1_ = false;
  abort_status2_ = false;

  // variance, projected quantities
  double f, fvar, df, dfvar;
  double raw_f, raw_fvar, raw_df, raw_dfvar;

  raw_f = en;
  raw_df = projected_grad(grad, search_dir);
  raw_fvar = variance(en_err);
  raw_dfvar = projected_gradvar(grad_err, search_dir);

  // initial f & df
  if (num_steps_ == 0) {
    f0_ = raw_f;
    df0_ = std::abs(raw_df);
  }

  // scaled obs
  std::tie(f,df,fvar,dfvar) = scale_obs(raw_f,raw_df,raw_fvar,raw_dfvar);
  /*
  std::cout << "\nstep = " << num_steps_ << "\n";
  std::cout << "f_raw, f = " << raw_f << ", "<< f  << "\n";
  std::cout << "df_raw, df = " << raw_df << ", "<< df  << "\n\n";
  */

  // Add the first observation to the gp
  GP_.add_obs(0.0, f, df, fvar, dfvar);

  // initial 'dt' value is 1.0
  double dt = rescale_t(1.0);
  current_t_ = 1.0;
  new_search_ = true;

  return dt;
}

double ProbLineSearch::do_step(const double& en, const double& en_err, 
  const RealVector& grad, const RealVector& grad_err, const RealVector& search_dir,
  bool& accept_status, std::string& message)
{
  if (num_steps_==0) {
    if (!new_search_) {
      throw std::range_error("ProbLineSearch::do_step: `start' method not called");
    }
    else new_search_ = false;
  }
  else { 
    if (new_search_) {
      throw std::range_error("ProbLineSearch::do_step: new `start' needed after accepted step");
    }
  }
  num_steps_++;

  // variance, projected quantities
  double f, fvar, df, dfvar;
  double raw_f, raw_fvar, raw_df, raw_dfvar;
  raw_f = en;
  raw_df = projected_grad(grad, search_dir);
  raw_fvar = variance(en_err);
  raw_dfvar = projected_gradvar(grad_err, search_dir);

  // If need to be aborted due to bad values
  if (std::isnan(raw_f) || std::isinf(raw_f) || 
     std::isnan(raw_df) || std::isinf(raw_df)) 
  {
    raw_f = 100.0;
    raw_df = 10.0;
    abort_status1_ = true;
  } 

  // scaled obs
  std::tie(f,df,fvar,dfvar) = scale_obs(raw_f,raw_df,raw_fvar,raw_dfvar);

  /*
  std::cout << "step = " << num_steps_ << "\n";
  std::cout << "f_raw, f = " << raw_f << ", "<< f  << "\n";
  std::cout << "df_raw, df = " << raw_df << ", "<< df  << "\n";
  getchar();
  */

  // Add the observation to the gp
  GP_.add_obs(current_t_, f, df, fvar, dfvar);
  GP_.update();

  // find next 't'
  double next_t, dt;
  if (abort_status1_ || num_steps_ >= max_steps_ ||
       num_expl_ >= max_expl_ || GP_.dmu(0.0) >= max_dmu0_) 
  {
    next_t = find_abort_t();
    abort_status2_ = true;
  }
  else {
    next_t = find_next_t();
  }

  dt = rescale_t(next_t-current_t_);
  current_t_ = next_t;

  accept_status = check_acceptibility();
  if (accept_status == true) {
    if (abort_status2_) message = "aborted";
    else message = "ok";

    /*
    # If this accept was not due to an abort and the step size did not change
    # *too much*, we use the accepted alpha as the new base step size alpha0
    # (and update a running average alpha_stats). Otherwise, we use said
    # running average as the new base step size.
    */
    /*
    double mcf = max_change_factor_; 
    double alpha = rescale_t(current_t_);
    bool alpha_good = false;
    if (alpha0_/mcf<alpha && alpha<alpha0_*mcf) alpha_good = true;
    if (!abort_status2_ && alpha_good) {
      alpha_stats_ = 0.95*alpha_stats_ + 0.05*alpha;
      alpha0_ = fpush_*alpha;
    }
    else {
      alpha0_ = alpha_stats_;
    }
    std::cout << "\n alpha0_ = " << alpha0_ << "\n";
    */

    // start new search next
    new_search_ = true; 
  }
  else {
    message = "continue";
  }

  return dt;
}

bool ProbLineSearch::check_acceptibility(void) 
{
  if (num_steps_ == 0) return false;

  if (abort_status2_) return true; 

  // Check Wolfe probability
  double pw = compute_p_wolfe(current_t_);
  if (pw > cW_) return true;
  else return false;
}

double ProbLineSearch::find_abort_t(void) 
{
  return 0.01;
}

double ProbLineSearch::find_next_t(void)
{
  // Find the step size for the next evaluation."""
  assert(num_steps_ >= 1);

  // Generate candidates: the points where the derivative of the posterior
  // mean equals the target value plus one exploration point to the right.
  //std::cout << " >>>>> in find_next_t\n";
  auto candidates = GP_.find_dmu_equal(target_df_);
  if (expl_policy_=="linear") {
    candidates.push_back(2.0*(num_expl_+1));
  }
  else if (expl_policy_=="exponential") {
    candidates.push_back(std::pow(2.0,(num_expl_+1)));
  }
  else {
    throw std::range_error("ProbLineSearch::find_next_t: unknown exploration policy");
  }
  // Compute p_Wolfe for candidates
  std::vector<double> pws;
  //int i = 0;
  //std::cout<<"\n"; 
  for (const double& t : candidates) {
    double pw = compute_p_wolfe(t);
    pws.push_back(pw);
    //std::cout<<i++<<" "<<t<<" "<<pw<<"\n"; 
  }
  int ibest = 0;
  double pwmax = pws[0]; 
  for (int i=1; i<pws.size(); ++i) {
    if (pwmax < pws[i]) {
      pwmax = pws[i];
      ibest = i;
    }
  }
  //std::cout<<"ibest = "<<ibest<<"\n"; 
  //getchar();

  // memorize when we have chosen the exploration point
  if (ibest == pws.size()-1) {
    num_expl_++;
  } 

  // Return the candidate t with maximal utility
  return candidates[ibest];
}

double ProbLineSearch::compute_p_wolfe(const double& t)
{
  // omputes the probability that step size ``t`` satisfies the adjusted
  //  Wolfe conditions under the current GP model.

  // Compute mean and covariance matrix of the two Wolfe quantities a and b
  // (equations (11) to (13) in [1]).
  double mu0 = GP_.mu(0.0);
  double dmu0 = GP_.dmu(0.0);
  double mu = GP_.mu(t);
  double dmu = GP_.dmu(t);
  double V0 = GP_.V(0.0);
  double Vd0 = GP_.Vd(0.0);
  double dVd0 = GP_.dVd(0.0);
  double dCov0t = GP_.dCov_0(t);
  double Covd0t = GP_.Covd_0(t);

  double x = c1_*t;
  double ma = mu0 - mu + x*dmu0;
  double Vaa = V0 + dVd0*x*x+GP_.V(t)+2.0*c1_*t*(Vd0-dCov0t)-2.0*GP_.Cov_0(t);
  double mb = dmu;
  double Vbb = GP_.dVd(t);

  // Very small variances can cause numerical problems. Safeguard against
  // this with a deterministic evaluation of the Wolfe conditions.
  if (Vaa<1.0E-9 || Vbb<1.0E-9) {
    if (ma>=0.0 && mb>=0.0) return 1.0;
    else return 0.0;
  }

  double Vab = Covd0t + x*GP_.dCovd_0(t) - GP_.Vd(t);
  // Compute correlation factor and integration bounds for adjusted p_Wolfe
  // and return the result of the bivariate normal integral.
  double rho = Vab/std::sqrt(Vaa*Vbb);
  double al = -ma/std::sqrt(Vaa);
  double bl = (df_lo_ - mb)/std::sqrt(Vbb);
  double bu = (df_hi_ - mb)/std::sqrt(Vbb);

  double inf = std::numeric_limits<double>::infinity();
  return bounded_bivariate_normal_integral(rho,al,inf,bl,bu);
}

double ProbLineSearch::bounded_bivariate_normal_integral(const double& rho, const double& xl, 
  const double& xu, const double& yl, const double& yu) const
{
  /*
  Computes the bounded bivariate normal integral.

  Computes the probability that ``xu >= X >= xl and yu >= Y >= yl`` where X
  and Y are jointly Gaussian random variables, with mean ``[0., 0.]`` and
  covariance matrix ``[[1., rho], [rho, 1.]]``.

  Inputs:
      :rho: Correlation coefficient of the bivariate normal random variable
      :xl, yl: Lower bounds of the integral
      :xu, yu: Upper bounds of the integral
  */
  double x1 = unbounded_bivariate_normal_integral(rho,xl,yl);
  double x2 = unbounded_bivariate_normal_integral(rho,xu,yl);
  double x3 = unbounded_bivariate_normal_integral(rho,xl,yu);
  double x4 = unbounded_bivariate_normal_integral(rho,xu,yu);
  double p = x1-x2-x3+x4;
  return std::max(0.0, std::min(p,1.0));
}

double ProbLineSearch::unbounded_bivariate_normal_integral(const double& rho_inp, 
  const double& xl, const double& yl) const
{
  /*
   Computes the unbounded bivariate normal integral.

   Computes the probability that ``X>=xl and Y>=yl`` where X and Y are jointly
    Gaussian random variables, with mean ``[0., 0.]`` and covariance matrix
  ``[[1., rho], [rho, 1.]]``.
  */
  double rho = std::max(-1.0,std::min(1.0,rho_inp));
  if (std::isinf(xl) || std::isinf(yl)) {
    return 0.0;
  }
  else if (std::isinf(-xl)) {
    if (std::isinf(-yl)) return 1.0;
    else CDF(-yl);
  }
  else if (std::isinf(-yl)) {
    return CDF(-xl);
  }
  else if (std::abs(rho)<1.0E-9) {
    return CDF(-xl)*CDF(-yl);
  }

  double PI = 3.1415926535897932384626433832795028841971693993751058209;
  double tp = 2.0*PI;
  double h = xl;
  double k = yl;
  double hk = h*k;
  double bvn = 0.0;

  // Gauss Legendre points and weights
  std::vector<double> ww;
  std::vector<double> xx;
  if (std::abs(rho) < 0.3) {
    // n = 6
    ww = {0.1713244923791705, 0.3607615730481384, 0.4679139345726904};
    xx = {0.9324695142031522, 0.6612093864662647, 0.2386191860831970};
  }
  else if (std::abs(rho) < 0.75) {
    // n = 12
    ww = {0.04717533638651177, 0.1069393259953183, 0.1600783285433464,
         0.2031674267230659, 0.2334925365383547, 0.2491470458134029};
    xx = {0.9815606342467191, 0.9041172563704750, 0.7699026741943050,
         0.5873179542866171, 0.3678314989981802, 0.1252334085114692};
  }
  else {
    // n = 20;
    ww = {0.01761400713915212, 0.04060142980038694, 0.06267204833410906,
         0.08327674157670475, 0.1019301198172404,  0.1181945319615184,
         0.1316886384491766,  0.1420961093183821,  0.1491729864726037,
         0.1527533871307259};
    xx = {0.9931285991850949, 0.9639719272779138, 0.9122344282513259,
         0.8391169718222188, 0.7463319064601508, 0.6360536807265150,
         0.5108670019508271, 0.3737060887154196, 0.2277858511416451,
         0.07652652113349733};
  }
  // full set
  std::vector<double> w;
  std::vector<double> x;
  for (int i=0; i<ww.size(); ++i) {
    w.push_back(ww[i]);
    x.push_back(1.0-xx[i]);
  }
  for (int i=0; i<ww.size(); ++i) {
    w.push_back(ww[i]);
    x.push_back(1.0+xx[i]);
  }

  if (std::abs(rho) < 0.925) {
    double hs = 0.5*(h*h + k*k);
    double asr = 0.5*std::asin(rho);
    std::vector<double> fx;
    for (const auto& xv : x) {
      double sn = std::sin(asr*xv);
      double sn_sq = sn*sn;
      double val = std::exp((sn*hk-hs)/(1.0-sn_sq));
      fx.push_back(val);
    }
    bvn = 0.0;
    for (int i=0; i<w.size(); ++i) {
      bvn += w[i]*fx[i];
    }
    bvn = bvn*asr/tp + CDF(-h)*CDF(-k);
  }
  else {
    if (rho < 0.0) {
      k = -k;
      hk = -hk;
    }
    if (std::abs(rho) < 1.0) {
      double ass = 1.0-rho*rho;
      double a = std::sqrt(ass);
      double bs = (h-k)*(h-k); 
      double asr = -0.5*(bs/ass + hk);
      double c = (4.0-hk)/8.0;
      double d = (12.0-hk)/80.0; 
      if (asr > -100.0) {
        bvn = a*std::exp(asr)*(1.0-c*(bs-ass)*(1.0-d*bs)/3.0 + c*d*ass*ass);
      }
      if (hk > -100.0) {
        double b = std::sqrt(bs);
        double sp = std::sqrt(tp)*CDF(-b/a);
        bvn = bvn - std::exp(-0.5*hk)*sp*b*(1.0-c*bs*(1.0-d*bs)/3.0);
      }
      a = 0.5*a;
      std::vector<double> xs;
      for (const auto& xv : x) {
        double ax = a*xv;
        xs.push_back(ax*ax);
      }
      double sum = 0.0;
      for (int i=0; i<x.size(); ++i) {
        double ax = a*x[i];
        double xs = ax*ax;
        double asr = -0.5*(bs/xs + hk);
        if (asr > -100.0) {
          double sp = 1.0 + c*xs*(1.+5.*d*xs);
          double rs = std::sqrt(1.0-xs);
          double rs2 = (1.0+rs)*(1.0+rs);
          double ep = std::exp(-0.5*hk*xs/rs2)/rs;
          double term = (a*std::exp(asr)*(sp-ep)*w[i]-bvn)/tp;
          sum += term;
        }
      }
      bvn = sum;
    }
    if (rho > 0.0) {
      bvn += CDF(-std::max(h, k));
    }
    else if (h >= k) {
      bvn = -bvn;
    }
    else {
      double L;
      if (h < 0.0) {
        L = CDF(k)-CDF(h);
      }
      else {
        L = CDF(-h)-CDF(-k);
      }
      bvn = L-bvn;
    }
  }

  return std::max(0.0, std::min(1.0, bvn));
}

double ProbLineSearch::CDF(const double& x) const
{
  // Cumulative density function (CDF) of the standard normal distribution.
  return 0.5 * (1.0 + std::erf(x/std::sqrt(2.0)));
}

ProbLineSearch::fquad ProbLineSearch::scale_obs(const double& raw_f, const double& raw_df, 
    const double& raw_fvar, const double& raw_dfvar)
{
  double beta = df0_*alpha0_;
  double f = (raw_f-f0_)/beta;
  double df = raw_df/df0_;
  double fvar = raw_fvar/(beta*beta);
  double dfvar = raw_dfvar/(df0_*df0_);
  return std::make_tuple(f,df,fvar,dfvar);
}

double ProbLineSearch::rescale_t(const double& raw_t) const 
{
  return raw_t * alpha0_;
}


double ProbLineSearch::variance(const double& ferr)
{
  return ferr*ferr;
}

double ProbLineSearch::projected_grad(const RealVector& grad, const RealVector& search_dir)
{
  return grad.dot(search_dir);
}


double ProbLineSearch::projected_gradvar(const RealVector& grad_err, const RealVector& search_dir)
{
  double sum = 0.0;
  for (int i=0; i<grad_err.size(); ++i) {
    double x = grad_err(i)*search_dir(i);
    sum += x*x;
  }
  return sum; // variance, not stddev
}

/*
Gaussian process functionality needed for the probabilistic
line search algorithm.
*/

ProbLineSearchGP::ProbLineSearchGP(void)
{
  theta_ = 1.0;
  theta_sq_ = theta_*theta_;
  offset_ = 10.0; // 10.0
  N_ = 0;
  ts_.clear();
  fs_.clear();
  dfs_.clear();
  fvars_.clear();
  dfvars_.clear();
  ready_ = false;
  have_min_obs_ = false;
}

int ProbLineSearchGP::set_parameters(const double& theta, const double& offset)
{
  theta_ = theta;
  offset_ = offset;
  return 0;
}

int ProbLineSearchGP::reset(void)
{
  N_ = 0;
  ts_.clear();
  fs_.clear();
  dfs_.clear();
  fvars_.clear();
  dfvars_.clear();

  ready_ = false;
  have_min_obs_ = false;

  return 0;
}

int ProbLineSearchGP::add_obs(const double& t, const double& f, const double& df,
    const double& fvar, const double& dfvar)
{
  /*
    Add a new observation (t, f, df, simga2_f, sigma2_df) to the GP.
    This stores the observation internally, but does NOT yet set up and invert
    the Gram matrix. Add observations with repeated calls to this method, then
    call ``gp.update()`` to set up and invert the Gram matrix. Only then you 
    can perform inference (calls to ``gp.mu(t)``, ``gp.V(t)``, etc...).
  */

  N_++;
  ts_.push_back(t);
  fs_.push_back(f);
  dfs_.push_back(df);
  fvars_.push_back(fvar);
  dfvars_.push_back(dfvar);

  ready_ = false;
  have_min_obs_ = false;

  return 0;
}

int ProbLineSearchGP::update(void)
{
  /*
  Set up the Gram matrix and compute its LU decomposition to make the GP
  ready for inference (calls to ``.gp.mu(t)``, ``gp.V(t)``, etc...).
    
  Call this method after you have manipulated the GP by
     - ``gp.reset()`` ing,
     - adding observations with ``gp.add(t, f, df)``, or
     - adjusting the sigmas via ``gp.update_sigmas()``.
  and want to perform inference next."""
  */

  if (ready_) return 0;

  // Set up the kernel matrices
  kernel_K_.resize(N_,N_);
  kernel_Kd_.resize(N_,N_);
  kernel_dKd_.resize(N_,N_);

  for (int i=0; i<N_; ++i) {
    for (int j=0; j<N_; ++j) {
      kernel_K_(i,j) = get_k(ts_[i],ts_[j]);
      kernel_Kd_(i,j) = get_kd(ts_[i],ts_[j]);
      kernel_dKd_(i,j) = get_dkd(ts_[i],ts_[j]);
    }
  }

  // the Gram matrix
  G_.resize(2*N_,2*N_);
  G_.block(0,0,N_,N_) = kernel_K_;
  G_.block(0,N_,N_,N_) = kernel_Kd_;
  G_.block(N_,0,N_,N_) = kernel_Kd_.transpose();
  G_.block(N_,N_,N_,N_) = kernel_dKd_;
  for (int i=0; i<N_; ++i) G_(i,i) += fvars_[i];
  for (int i=0; i<N_; ++i) G_(N_+i,N_+i) += dfvars_[i];

  // LU decomposition of G
  //FullPivLU_.compute(G_);

  // Pre-compute the regression weights used in mu
  w_.resize(2*N_);
  RealVector b(2*N_);
  for (int i=0; i<N_; ++i) b(i) = fs_[i];
  for (int i=0; i<N_; ++i) b(N_+i) = dfs_[i];
  //w_ = FullPivLU_.solve(b);
  w_ = G_.fullPivLu().solve(b);

  // ready now
  ready_ = true;

  return 0;
}

RealVector ProbLineSearchGP::solve_G(const RealVector& b) const
{
  assert(ready_);
  //return FullPivLU_.solve(b); 
  return G_.fullPivLu().solve(b);
}

double ProbLineSearchGP::mu(const double& t) const
{
  assert(ready_);
  // Evaluate posterior mean of f at ``t``
  RealVector kvec(2*N_);
  for (int i=0; i<N_; ++i) kvec[i] = get_k(t, ts_[i]);
  for (int i=0; i<N_; ++i) kvec[N_+i] = get_kd(t, ts_[i]);
  return w_.dot(kvec);
}

double ProbLineSearchGP::dmu(const double& t) const
{
  assert(ready_);
  // Evalulate first derivative of the posterior mean of df at ``t``
  RealVector kvec(2*N_);
  for (int i=0; i<N_; ++i) kvec[i] = get_kd(ts_[i],t);
  for (int i=0; i<N_; ++i) kvec[N_+i] = get_dkd(t, ts_[i]);
  return w_.dot(kvec);
}

double ProbLineSearchGP::d2mu(const double& t) const
{
  assert(ready_);
  // Evaluate 2nd derivative of the posterior mean of f at ``t``
  RealVector kvec(2*N_);
  for (int i=0; i<N_; ++i) kvec[i] = get_d2k(t,ts_[i]);
  for (int i=0; i<N_; ++i) kvec[N_+i] = get_d2kd(t,ts_[i]);
  return w_.dot(kvec);
}

double ProbLineSearchGP::d3mu(const double& t) const
{
  assert(ready_);
  // Evaluate 3rd derivative of the posterior mean of f at ``t``
  RealVector kvec(2*N_);
  for (int i=0; i<N_; ++i) kvec[i] = get_d3k(t,ts_[i]);
  for (int i=0; i<N_; ++i) kvec[N_+i] = 0.0;
  return w_.dot(kvec);
}

double ProbLineSearchGP::V(const double& t) const
{
  // Evaluate posterior variance of f at ``t``
  assert(ready_);
  RealVector kvec(2*N_);
  for (int i=0; i<N_; ++i) kvec[i] = get_k(t,ts_[i]);
  for (int i=0; i<N_; ++i) kvec[N_+i] = get_kd(t,ts_[i]);
  double ktt = get_k(t,t);

  return ktt - kvec.dot(solve_G(kvec));
}

double ProbLineSearchGP::Vd(const double& t) const
{
  // Evaluate posterior co-variance of f and df at ``t``
  assert(ready_);
  RealVector kvec_a(2*N_);
  for (int i=0; i<N_; ++i) kvec_a[i] = get_k(t,ts_[i]);
  for (int i=0; i<N_; ++i) kvec_a[N_+i] = get_kd(t,ts_[i]);
  RealVector kvec_b(2*N_);
  for (int i=0; i<N_; ++i) kvec_b[i] = get_kd(ts_[i],t);
  for (int i=0; i<N_; ++i) kvec_b[N_+i] = get_dkd(t,ts_[i]);
  double ktt = get_k(t,t);

  return ktt - kvec_a.dot(solve_G(kvec_b));
}

double ProbLineSearchGP::dVd(const double& t) const
{
  // Evaluate posterior variance of df at ``t`
  assert(ready_);
  RealVector kvec(2*N_);
  for (int i=0; i<N_; ++i) kvec[i] = get_kd(ts_[i],t);
  for (int i=0; i<N_; ++i) kvec[N_+i] = get_dkd(t,ts_[i]);
  double dkdtt = get_dkd(t,t);

  return dkdtt - kvec.dot(solve_G(kvec));
}

double ProbLineSearchGP::Cov_0(const double& t) const
{
  // Evaluate posterior co-variance of f at 0. and ``t``
  assert(ready_);
  RealVector kvec_a(2*N_);
  for (int i=0; i<N_; ++i) kvec_a[i] = get_k(0.0,ts_[i]);
  for (int i=0; i<N_; ++i) kvec_a[N_+i] = get_kd(0.0,ts_[i]);
  RealVector kvec_b(2*N_);
  for (int i=0; i<N_; ++i) kvec_b[i] = get_k(t,ts_[i]);
  for (int i=0; i<N_; ++i) kvec_b[N_+i] = get_kd(t,ts_[i]);
  double k0t = get_k(0.0,t);

  return k0t - kvec_a.dot(solve_G(kvec_b));
}

double ProbLineSearchGP::Covd_0(const double& t) const
{
  // Evaluate posterior co-variance of f at 0. and df at ``t``
  // !!! I changed this in line_search new, Covd_0 <-> dCov_0
  assert(ready_);
  RealVector kvec_a(2*N_);
  for (int i=0; i<N_; ++i) kvec_a[i] = get_k(0.0,ts_[i]);
  for (int i=0; i<N_; ++i) kvec_a[N_+i] = get_kd(0.0,ts_[i]);
  RealVector kvec_b(2*N_);
  for (int i=0; i<N_; ++i) kvec_b[i] = get_kd(ts_[i],t);
  for (int i=0; i<N_; ++i) kvec_b[N_+i] = get_dkd(t,ts_[i]);
  double kd0t = get_kd(0.0,t);

  return kd0t - kvec_a.dot(solve_G(kvec_b));
}

double ProbLineSearchGP::dCov_0(const double& t) const
{
  // Evaluate posterior co-variance of df at 0. and f at ``t``
  // !!! I changed this in line_search new, Covd_0 <-> dCov_0
  assert(ready_);
  RealVector kvec_a(2*N_);
  for (int i=0; i<N_; ++i) kvec_a[i] = get_kd(ts_[i],0.0);
  for (int i=0; i<N_; ++i) kvec_a[N_+i] = get_dkd(0.0,ts_[i]);
  RealVector kvec_b(2*N_);
  for (int i=0; i<N_; ++i) kvec_b[i] = get_k(t,ts_[i]);
  for (int i=0; i<N_; ++i) kvec_b[N_+i] = get_kd(t,ts_[i]);
  double dk0t = get_kd(t,0.0);

  return dk0t - kvec_a.dot(solve_G(kvec_b));
}

double ProbLineSearchGP::dCovd_0(const double& t) const
{
  // Evaluate posterior co-variance of df at 0. and ``t`
  assert(ready_);
  RealVector kvec_a(2*N_);
  for (int i=0; i<N_; ++i) kvec_a[i] = get_kd(ts_[i],0.0);
  for (int i=0; i<N_; ++i) kvec_a[N_+i] = get_dkd(0.0,ts_[i]);
  RealVector kvec_b(2*N_);
  for (int i=0; i<N_; ++i) kvec_b[i] = get_kd(ts_[i],t);
  for (int i=0; i<N_; ++i) kvec_b[N_+i] = get_dkd(t,ts_[i]);
  double dkd0t = get_dkd(0.0,t);

  return dkd0t - kvec_a.dot(solve_G(kvec_b));
}

int ProbLineSearchGP::cubic_poly_coeff(const double& t, 
  double& a, double& b, double& c, double& d) const
{
  // The posterior mean ``mu`` of this GP is piece-wise cubic. Return the
  //  coefficients of the cubic polynomial that is ``mu`` at ``t``."""
  double d1 = dmu(t);
  double d2 = d2mu(t);
  double d3 = d3mu(t);
  a = d3/6.0;
  b = 0.5*d2-3*a*t;
  c = d1-3*a*t*t-2*b*t;
  d = mu(t)-a*t*t*t-b*t*t-c*t;

  return 0;
}

int ProbLineSearchGP::quad_poly_coeff(const double& t, double& a, double& b, double& c) const
{
  /* 
    The posterior mean ``mu`` of this GP is piece-wise cubic. Return the
    coefficients of the **quadratic** polynomial that is the **derivative** of
    ``mu`` at ``t`
  */
  double d1 = dmu(t);
  double d2 = d2mu(t);
  double d3 = d3mu(t);
  a = 0.5*d3;
  b = d2 - d3*t;
  c = d1 - d2*t + 0.5*d3*t*t;

  return 0;
}

std::vector<double> ProbLineSearchGP::find_dmu_equal(const double& val) const
{
  // Finds points where the derivative of the posterior mean equals ``val``
  //   and the second derivative is positive.

  // We want to go through the observations from smallest to largest t
  std::vector<double> ts_sorted = ts_;
  std::sort(ts_sorted.begin(),ts_sorted.end());
  //std::cout << "t values\n";
  //for (const auto& t : ts_sorted) {
  //  std::cout << t << "\n";
  //}

  std::vector<double> solutions;
  std::vector<double> solutions_cell;
  for (int i=0; i<N_-1; ++i) {
    double t1 = ts_sorted[i];
    double t2 = ts_sorted[i+1];
    double t = t1+0.5*(t2-t1);
    double a,b,c;
    quad_poly_coeff(t,a,b,c);
    solutions_cell = quad_poly_solve(a,b,c,val);
    for (const double& r : solutions_cell) {
      if (r>t1 && r<t2) solutions.push_back(r);
    }
  }
  /*
  std::cout << "solutions\n";
  for (const auto& s : solutions) {
    std::cout << s << "\n";
  }
  //getchar();
  */

  return solutions;
}

std::vector<double> ProbLineSearchGP::find_cubic_minima(void) const
{
  /* Find the local minimizers of the posterior mean.

    The posterior mean is a  cubic polynomial in each of the cells"
    [t_i, t_i+1] where the t_i are the sorted observed ts. For each of these
    cells, return the minimizer of the cubic polynomial if it exists and
    happens to lie in that cell.
  */

  return find_dmu_equal(0.0);
}

double ProbLineSearchGP::expected_improvement(const double& t) 
{
  /* Computes the expected improvement at position ``t`` 
     under the current GP model. 
     Reference "current best" is the observed ``t`` with minimal 
     posterior mean.
  */

  // Find the observation with minimal posterior mean, if it has not yet been
  // computed by a previous call to this method
  if (!have_min_obs_) {
    std::vector<double> mean_obs;
    for (const double& t : ts_) mean_obs.push_back(mu(t));
    std::sort(mean_obs.begin(),mean_obs.end());
    min_posterior_mean_ = mean_obs[0];
    have_min_obs_ = true;
  }

  // Compute posterior mean and variance at t
  double m = mu(t);
  double v = V(t);
  double d = min_posterior_mean_-m;
  double t1 = 0.5 * d * (1.0+std::erf(d/std::sqrt(2.0*v)));
  double t2 = std::sqrt(0.5*v/PI_) * std::exp(-0.5*d*d/v);
  
  return t1+t2;  
}

std::vector<double> ProbLineSearchGP::quad_poly_solve(const double& a, 
  const double& b, const double& c, const double& rhs) const
{
  // Computes *real* solutions of f'(t) = a*t**2 + b*t + c = rhs with f''(t)>0
  //std::cout << "a, b, c, rhs ="<<a<<", "<<b<<", "<<c<<" "<<rhs<<"\n";
  //std::cout << "det ="<<b*b - 4.0*a*(c-rhs)<<"\n";

  std::vector<double> solutions;
  // Check if a is almost zero. If so, solve the remaining linear equation. Note
  // that we return only soultions with f''(t) = b > 0
  if (std::abs(a)<1.0E-9) {
    if (b > 1.0E-9) {
      solutions.push_back((rhs-c)/b);
    }
    return solutions;
  }

  // Compute the term under the square root in pq formula, if it is negative,
  // there is no real solution
  double det = b*b - 4.0*a*(c-rhs);
  if (det < 0.0) return solutions;

  // Otherwise, compute the two roots 
  //int sign_a;
  //if (a>=0.0) sign_a = 1;
  //else sign_a = -1;
  double s = std::sqrt(det);
  double r1 = (-b - s)/(2.0*a);
  double r2 = (-b + s)/(2.0*a);
  //std::cout << "rhs = "<<rhs<<"\n";
  //std::cout << "r1, r2 = "<<r1<<",  "<<r2<<"\n";

  // Return the one with f''(t) = 2at + b > 0, or []
  if (2*a*r1+b > 0.0) solutions.push_back(r1);
  else if (2*a*r2+b > 0.0) solutions.push_back(r2);

  return solutions;
}

double ProbLineSearchGP::get_k(const double& x, const double& y) const
{
  // kernel function
  double min = offset_ + std::min(x,y);
  double min_sq = min*min;
  return theta_sq_*(min_sq*min/3.0 + 0.5*std::abs(x-y)*min_sq);
}

double ProbLineSearchGP::get_kd(const double& x, const double& y) const
{
  // derivative of the kernel function, 1st derivative w.r.t. right argument
  double xx = x + offset_;
  double yy = y + offset_;
  if (x < y) return 0.5*theta_sq_*xx*xx;
  else return theta_sq_*(xx*yy - 0.5*yy*yy);
}

double ProbLineSearchGP::get_dkd(const double& x, const double& y) const
{
  // derivative of the kernel function, 1st derivative w.r.t. both arguments.
  double xx = x + offset_;
  double yy = y + offset_;
  return theta_sq_*std::min(xx,yy);
}

double ProbLineSearchGP::get_d2k(const double& x, const double& y) const
{
  // Derivative of kernel function,  2nd derivative w.r.t. left argument.
  if (x < y) return theta_sq_*(y-x);
  else return 0.0;
}

double ProbLineSearchGP::get_d3k(const double& x, const double& y) const
{
  // Derivative of kernel function,  3rd derivative w.r.t. left argument.
  if (x < y) return -theta_sq_;
  else return 0.0;
}

double ProbLineSearchGP::get_d2kd(const double& x, const double& y) const
{
  // Derivative of kernel function, 2nd derivative w.r.t. left argument,
  //  1st derivative w.r.t. right argument.
  if (x < y) return theta_sq_;
  else return 0.0;
}

int ProbLineSearchGP::test(void) 
{
  // quadratic polynomial solve
  auto sol = quad_poly_solve(1., 3., -2., -4.);
  std::cout << "sol : ["; 
  for (const auto& s : sol) std::cout << s << " ";
  std::cout << "] =? " << "[-1]\n"; 

  std::vector<double> ts;
  std::vector<double> fs;
  std::vector<double> dfs;
  std::default_random_engine gen;
  std::normal_distribution<double> dist;
  for (int i=0; i<10; ++i) {
    ts.push_back(dist(gen));
    fs.push_back(dist(gen));
    dfs.push_back(dist(gen));
  }
  for (int i=0; i<10; ++i) add_obs(ts[i],fs[i],dfs[i]);
  update();
  for (int i=0; i<10; ++i) {
    double t = ts[i];
    double f = fs[i];
    double df = dfs[i];
    std::cout << "V = " << V(t) << " test " << V(t)-1.0E-9<<"\n";
    std::cout << "mu = " << mu(t) << " test " << mu(t)-f<<"\n";
    std::cout << "dmu = " << dmu(t) << " test " << dmu(t)-df<<"\n\n";
  }
  return 0;
}

} // end namespace opt














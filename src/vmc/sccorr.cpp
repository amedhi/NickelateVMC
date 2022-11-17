/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-10 21:35:10
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-11 00:08:00
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./sccorr.h"
#include "../utils/curve_fit.h"


namespace vmc {

void SC_Correlation::setup(const lattice::Lattice& lattice, const var::MF_Order::pairing_t& pairsymm)
{
  // SC pairing symmetry
  pair_symm_ = pairsymm;

  // max distance (catch sites along a1)
  Vector3d R = lattice.size1() * lattice.vector_a1();
  max_dist_ = std::nearbyint(R.norm());
  max_dist_ = max_dist_/2;
  //std::cout << "max_dist = "<<max_dist_<<"\n";
  sitepair_list_.resize(max_dist_+1);

  for (int i=0; i<lattice.num_sites(); ++i) {
    Vector3i n1 = lattice.site(i).bravindex();
    Vector3d r1 = lattice.site(i).coord();
    for (int j=i; j<lattice.num_sites(); ++j) {
      Vector3i n2 = lattice.site(j).bravindex();
      // pick-up pairs of sites in unit cells along '+x'
      if (n1(1)!=n2(1) || n1(2)!=n2(2)) break;

      //  pairs of sites should lie in straight line along '+x'
      Vector3d R = r1-lattice.site(j).coord();
      if (std::abs(R(1))>1.0E-6 or std::abs(R(2))>1.0E-6) continue;

      // select pairs at ditance <= max_dist_
      int d = std::nearbyint(std::abs(R(0)));
      if (d > max_dist_) continue;
      sitepair_list_.pairs_at_dist(d).push_back({i,j});
      //std::cout<<i<<" -- "<<j<<" = "<<d<<"\n"; //getchar();
    }
  }

  // bond pair types for this different lattces 
  if (lattice.id()==lattice::lattice_id::NICKELATE) {
    bondpair_types_.resize(4);
    bondpair_types_[0] = std::make_pair(0,0);
    bondpair_types_[1] = std::make_pair(0,1);
    bondpair_types_[2] = std::make_pair(5,5);
    bondpair_types_[3] = std::make_pair(5,6);
  }
  else if (lattice.id()==lattice::lattice_id::NICKELATE_2D ||
           lattice.id()==lattice::lattice_id::NICKELATE_2L) {
    bondpair_types_.resize(4);
    bondpair_types_[0] = std::make_pair(0,0);
    bondpair_types_[1] = std::make_pair(0,1);
    bondpair_types_[2] = std::make_pair(3,3);
    bondpair_types_[3] = std::make_pair(3,4);
  }
  else {
    bondpair_types_.resize(2);
    bondpair_types_[0] = std::make_pair(0,0);
    bondpair_types_[1] = std::make_pair(0,1);
  }

  if (pair_symm_==var::MF_Order::pairing_t::SWAVE) {
    throw std::range_error("SC_Correlation::setup: not implemented for this symmetry");
    /*
    bondpair_types_.resize(num_basis_sites_);
    for (int i=0; i<num_basis_sites_; ++i) {
      bondpair_types_[i] = std::make_pair(i,i);
    }
    */
  }

  // allocate storages
  //corr_data_.resize(max_dist_+1, bondpair_types_.size());
  //count_.resize(max_dist_+1, bondpair_types_.size());

  // At location 'max_dist_+1': to store the value of 'phi'
  corr_data_.resize(max_dist_+2, bondpair_types_.size());
  count_.resize(max_dist_+2, bondpair_types_.size());

  std::vector<std::string> elem_names;
  for (const auto& p : bondpair_types_)
    elem_names.push_back(std::to_string(p.first)+"-"+std::to_string(p.second));
  this->resize(corr_data_.size(), elem_names);
}

void SC_Correlation::measure(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  if (pair_symm_==var::MF_Order::pairing_t::SWAVE) {
    measure_swave(lattice, model, config);
  }
  else if (pair_symm_==var::MF_Order::pairing_t::DWAVE) {
    measure_dwave(lattice, model, config);
  }
  else if (pair_symm_==var::MF_Order::pairing_t::EXTENDED_S) {
    measure_dwave(lattice, model, config);
  }
  else if (pair_symm_==var::MF_Order::pairing_t::null) {
    measure_dwave(lattice, model, config);
  }
  else {
    throw std::range_error("SC_Correlation::measure: not implemented for this symmetry");
  }
}

void SC_Correlation::measure_dwave(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  corr_data_.setZero();
  count_.setZero();
  double ampl;
  int fr_site_i, fr_site_ia, to_site_j, to_site_jb;
  for (int d=0; d<=max_dist_; ++d) {
    for (const auto& p : sitepair_list_.pairs_at_dist(d)) {
      // pairs of two sites
      fr_site_i  = p.first;
      to_site_j = p.second;

      for (const auto& id1 : lattice.site(p.first).outbond_ids()) {
        // skip boundary crossing bonds
        if (lattice.bond(id1).sign()<0) continue;

        int t1 = lattice.bond(id1).type();
        // check if SC pair defined for this bond 
        bool sc_pair = false;
        for (int m=0;  m<bondpair_types_.size(); ++m) {
          if (t1==bondpair_types_[m].first) {
            sc_pair = true; break;
          }
        }
        if (!sc_pair) continue;
        // bond-partner of first site
        fr_site_ia = lattice.bond(id1).tgt_id();

        // second bond
        for (const auto& id2 : lattice.site(p.second).outbond_ids()) {
          // skip boundary crossing bonds
          if (lattice.bond(id2).sign()<0) continue;

          int t2 = lattice.bond(id2).type();
          // bond-partner of second site
          to_site_jb = lattice.bond(id2).tgt_id();

          // calculate bond-singlet hopping amplitude
          for (int m=0;  m<bondpair_types_.size(); ++m) {
            if (t1==bondpair_types_[m].first && t2==bondpair_types_[m].second) {
              ampl = std::real(config.apply_bondsinglet_hop(fr_site_i,
                  fr_site_ia, to_site_j, to_site_jb));
              // Hermitian conjugate term
              ampl += std::real(config.apply_bondsinglet_hop(to_site_j,
                  to_site_jb, fr_site_i, fr_site_ia));
              corr_data_(d,m) += 0.5*ampl; // average of HC terms
              count_(d,m) += 1;
              //double term = std::real(config.apply_cdagc_up(i_cdag,j_c,1));
              //corr_data_(d,m) += term;
              //std::cout<<i_cdag<<"-"<<ia_cdag<<" x "<<j_c<<"-"<<ja_c<<"\n";
              //std::cout << "d = "<<d << " term = "<<term << "\n";
              //getchar();
            }
          }
        }
      }
    }
  }

  // average over number of bonds
  for (int d=0; d<=max_dist_; ++d) {
    for (int m=0;  m<bondpair_types_.size(); ++m) {
      if (count_(d,m) != 0) {
        corr_data_(d,m) /= count_(d,m);
      }
    }
  }

  // SC order parameters
  int d = max_dist_+1;
  for (int m=0;  m<bondpair_types_.size(); ++m) {
    double phi = std::sqrt(std::abs(corr_data_(max_dist_,m)));
    corr_data_(d,m) = phi;
  }

  // reshape add to databin
  *this << Eigen::Map<mcdata::data_t>(corr_data_.data(), corr_data_.size());
}

void SC_Correlation::measure_swave(const lattice::Lattice& lattice, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  corr_data_.setZero();
  int i_cdag, i_c;
  for (int d=1; d<max_dist_; ++d) {
    for (int n=0; n<num_basis_sites_; ++n) {
      for (const auto& p : symm_list_[n].pairs_at_dist(d)) {
        i_cdag = p.first;
        i_c = p.second;
        //std::cout << "swave corr ["<<n<<"] : "<<i_cdag<<"--"<<i_c<<"\n"; getchar();
        double term = std::real(config.apply_sitepair_hop(i_cdag,i_c));
        //corr_data_(d,n) += term;
        term += std::real(config.apply_sitepair_hop(i_c,i_cdag));
        corr_data_(d,n) += 0.5 * term;
      }
    }
  }
  for (int d=1; d<max_dist_; ++d) {
    for (int n=0; n<num_basis_sites_; ++n) {
      corr_data_(d,n) /= symm_list_[n].pairs_at_dist(d).size();
    }
  }

  //std::cout << max_dist_ << "\n\n";
  //std::cout << corr_data_ << "\n\n";
  //std::cout << Eigen::Map<mcdata::data_t>(corr_data_.data(), corr_data_.size()) << "\n";
  //getchar();

  // reshape add to databin
  *this << Eigen::Map<mcdata::data_t>(corr_data_.data(), corr_data_.size());
}


void SC_Correlation::print_heading(const std::string& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!is_open()) open_file();
  fs_ << header;
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  // total value
  for (const auto& p : xvars) fs_ << std::setw(14)<<p.substr(0,14);
  fs_ << std::setw(6)<<"d";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  /*if (xvars.size()>0) {
    fs_ << "# ";
    for (const auto& p : xvars) fs_<<std::setw(14)<<p.substr(0,14);
    fs_ << std::endl;
  }*/
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void SC_Correlation::print_result(const std::vector<double>& xvals) 
{
  if (!is_on()) return;
  if (!is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  /*if (xvals.size()>0) {
    fs_ << "#";
    for (const auto& p : xvals) fs_ << std::setw(14) << p;
    fs_ << std::endl;
  }*/

  // reshape back
  //corr_data_ = Eigen::Map<Eigen::MatrixXd>(MC_Data::mean_data(), corr_data_.size());

  std::vector<double> odlro;
  for (int d=0; d<=max_dist_; ++d) {
    for (const auto& p : xvals) fs_ << std::setw(14) << p;
    fs_ << std::setw(6) << d; 
    int n = d;
    for (int i=0; i<bondpair_types_.size(); ++i) {
      fs_ << MC_Data::result_str(n);
      n += (max_dist_+2);
    }
    fs_ << MC_Data::conv_str(d); 
    fs_ << std::endl; 
  }

  // Order parameter (from max_distance value)
  fs_ << "\n# phi(1)"; 
  for (const auto& p : xvals) fs_ << std::scientific << std::setw(14) << p;
  //fs_ << std::setw(6) << max_dist_-1; 
  int n=max_dist_+1;
  for (int i=0; i<bondpair_types_.size(); ++i) {
    fs_ << MC_Data::result_str(n);
    n += (max_dist_+2);
  }
  fs_ << MC_Data::conv_str(max_dist_+1); 
  fs_ << std::endl; 

  // Order parameter 'phi': obtained by fitting the data
  fs_ << "# phi(2)"; 
  for (const auto& p : xvals) fs_ << std::scientific << std::setw(14) << p;
  util::CurveFit curve_fit;
  for (int i=0; i<bondpair_types_.size(); ++i) {
    Eigen::VectorXd xv(max_dist_-2);
    Eigen::VectorXd yv(max_dist_-2);
    Eigen::VectorXd yerr(max_dist_-2);
    int j = 0;
    int n = max_dist_*i;
    for (int d=2; d<max_dist_; ++d) {
      xv(j) = d;
      yv(j) = MC_Data::mean(n+d);
      yerr(j) = MC_Data::stddev(n+d);
      j++;
    }
    //std::cout << yv.transpose() << "\n";
    //std::cout << yerr.transpose() << "\n\n";
    // fit
    double phi, err;
    Eigen::VectorXd p = Eigen::VectorXd::Ones(4);
    Eigen::VectorXd p2 = Eigen::VectorXd::Ones(3);
    FitFunc func;
    FitFunc2 func2;
    if (curve_fit.lmfit(func,p,xv,yv)) {
      phi = std::sqrt(std::abs(p(0)));
    }
    else if (curve_fit.lmfit(func2,p2,xv,yv)) {
      phi = std::sqrt(std::abs(p2(0)));
    }
    else {
      phi = std::sqrt(std::abs(yv(max_dist_-1)));
    }
    err = std::sqrt(yerr.sum()/yerr.size());
    fs_ << std::scientific << std::setw(15) << phi;
    fs_ << std::fixed << std::setw(10) << err;
  }
  fs_ << MC_Data::conv_str(0); 
  fs_ << std::endl; 

  //fstream() << MC_Data::conv_str(0); //.substr(0,10); 
  fs_ << std::flush;
  close_file();
}

} // end namespace vmc


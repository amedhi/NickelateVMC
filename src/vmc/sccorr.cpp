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

void SC_Correlation::setup(const lattice::LatticeGraph& graph)
{
  lattice::LatticeGraph::out_edge_iterator ei, ei_end;
  max_dist_ = graph.lattice().size1()/2+1;
  num_basis_sites_ = graph.lattice().num_basis_sites();
  // for each basis site, pairs of sites connected by translational symmetry
  //std::cout << "--SC_Correlation::setup: HACK--\n";
  symm_list_.resize(num_basis_sites_); 
  for (auto& elem : symm_list_) elem.resize(max_dist_);
  for (auto s1=graph.sites_begin(); s1!=graph.sites_end(); ++s1) {
  //for (auto s1=graph.sites_begin(); s1!=graph.sites_begin()+1; ++s1) {
    Vector3d rs1 = graph.site_cellcord(s1);
    //if (rs1[0]>0) continue;
    for (auto s2=s1; s2!=graph.sites_end(); ++s2) {
      if (graph.site_uid(s1)!=graph.site_uid(s2)) continue;
      Vector3d rs2 = graph.site_cellcord(s2);
      // pick sites along direction-1 (x-dir)
      if (rs1[1]!=rs2[1] ) break;
      if (rs1[2]!=rs2[2] ) break;
      int d = std::abs(rs1[0]-rs2[0]);
      int n = graph.site_uid(s2);
      //std::cout << graph.site(s1) << "  " << graph.site(s2) << "\n";
      if (d<max_dist_) {
        symm_list_[n].pairs_at_dist(d).push_back({graph.site(s1), graph.site(s2)});
      }
    }
  }

  /*
  for (int n=0; n<num_basis_sites_; ++n) {
    for (int d=1; d<max_dist_; ++d) {
      std::cout << "pairs = " << symm_list_[n].pairs_at_dist(d).size() << "\n";
      for (const auto& p : symm_list_[n].pairs_at_dist(d)) {
        std::cout << "(n,d) ="<<n<<", "<<d<<": ";
        std::cout << p.first<<"--"<<p.second<<"\n";
        getchar();
      }
    }
  }
  getchar();
  */

  // bond pair types for this lattice
  if (graph.lattice().id()==lattice::lattice_id::NICKELATE) {
    bondpair_types_.resize(4);
    bondpair_types_[0] = std::make_pair(0,0);
    bondpair_types_[1] = std::make_pair(0,1);
    bondpair_types_[2] = std::make_pair(5,5);
    bondpair_types_[3] = std::make_pair(5,6);
  }
  else if (graph.lattice().id()==lattice::lattice_id::NICKELATE_2D ||
           graph.lattice().id()==lattice::lattice_id::NICKELATE_2L) {
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
  //corr_data_.resize(bondpair_types_.size(), max_dist_);
  corr_data_.resize(max_dist_, bondpair_types_.size());
  std::vector<std::string> elem_names;
  for (const auto& p : bondpair_types_)
    elem_names.push_back(std::to_string(p.first)+"-"+std::to_string(p.second));
  this->resize(corr_data_.size(), elem_names);
}

void SC_Correlation::measure(const lattice::LatticeGraph& graph, 
  const model::Hamiltonian& model, const SysConfig& config)
{
  lattice::LatticeGraph::out_bond_iterator b1, b1_end, b2, b2_end;
  corr_data_.setZero();
  int i_cdag, ia_cdag, j_c, ja_c, i_phase, j_phase;
  for (int d=1; d<max_dist_; ++d) {
    for (int n=0; n<num_basis_sites_; ++n) {
      for (const auto& p : symm_list_[n].pairs_at_dist(d)) {
        for (std::tie(b1,b1_end)=graph.out_bonds(p.first); b1!=b1_end; ++b1) {
          int t1 = graph.bond_type(b1);
          for (std::tie(b2,b2_end)=graph.out_bonds(p.second); b2!=b2_end; ++b2) {
            int t2 = graph.bond_type(b2);
            for (int m=0;  m<bondpair_types_.size(); ++m) {
              if (t1==bondpair_types_[m].first && t2==bondpair_types_[m].second) {
                i_cdag = p.first;
                ia_cdag = graph.target(b1);
                i_phase = graph.bond_sign(b1);
                j_c = p.second;
                ja_c = graph.target(b2);
                j_phase = graph.bond_sign(b2);
                double term = std::real(config.apply_bondsinglet_hop(i_cdag,ia_cdag,
                              i_phase,j_c,ja_c,j_phase));
                corr_data_(d,m) += term;
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
  }
  // average over number of bonds
  if (graph.lattice().id()==lattice::lattice_id::NICKELATE ||
      graph.lattice().id()==lattice::lattice_id::NICKELATE_2D ||
      graph.lattice().id()==lattice::lattice_id::NICKELATE_2L) {
    for (int d=1; d<max_dist_; ++d) {
      corr_data_(d,0) /= symm_list_[0].pairs_at_dist(d).size();
      corr_data_(d,1) /= symm_list_[0].pairs_at_dist(d).size();
      corr_data_(d,2) /= symm_list_[1].pairs_at_dist(d).size();
      corr_data_(d,3) /= symm_list_[1].pairs_at_dist(d).size();
    }
  }
  else {
    for (int d=1; d<max_dist_; ++d) {
      for (int m=0;  m<bondpair_types_.size(); ++m) {
        corr_data_(d,m) /= symm_list_[0].pairs_at_dist(d).size();
      }
    }
  }
  //std::cout << corr_data_ << "\n"; getchar();

  // reshape add to databin
  *this << Eigen::Map<mcdata::data_t>(corr_data_.data(), corr_data_.size());

  /*
  for (int d=0; d<max_dist_; ++d) {
    bond_pair_corr_[d].setZero();
  }
  for (int i=0; i<src_pairs_.size(); ++i) {
    int d = pair_distance_[i];  
    for (std::tie(b1,b1_end)=graph.out_bonds(src_pairs_[i].first); b1!=b1_end; ++b1) {
      unsigned i_cdag = graph.source(b1);
      unsigned ia_cdag = graph.target(b1);
      unsigned type_i = graph.bond_type(b1);
      int i_phase = graph.bond_sign(b1);
      for (std::tie(b2,b2_end)=graph.out_bonds(src_pairs_[i].second); b2!=b2_end; ++b2) {
        unsigned j_c = graph.source(b2);
        unsigned ja_c = graph.target(b2);
        unsigned type_j = graph.bond_type(b2);
        int j_phase = graph.bond_sign(b2);
        double term = std::real(config.apply_bondsinglet_hop(i_cdag,ia_cdag,
          i_phase,j_c,ja_c,j_phase));
        // pair corr
        bond_pair_corr_[d](type_i,type_j) += term;
        //std::cout << "d = "<<d << " term = "<<term << "\n"; 
      }
    }
  }
  int n = max_dist_ * num_bond_types_ * num_bond_types_;
  mcdata::data_t corr_data_(n);
  n = 0;
  for (int d=0; d<max_dist_; ++d) {
    for (int i=0; i<num_bond_types_; ++i) {
      for (int j=0; j<num_bond_types_; ++j) {
        int multi = std::max(1,num_symm_pairs_[d](i,j));
        corr_data_[n] = bond_pair_corr_[d](i,j)/multi;
        ++n;
      }
    }
  }
  // add to databin
  *this << corr_data_;
  */
  //std::cout << "SC_Correlation::measure\n"; getchar();
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
  std::vector<double> odlro;
  for (int d=2; d<max_dist_; ++d) {
    for (const auto& p : xvals) fs_ << std::setw(14) << p;
    fs_ << std::setw(6) << d; 
    int n = d;
    for (int i=0; i<bondpair_types_.size(); ++i) {
      fs_ << MC_Data::result_str(n);
      if (d==(max_dist_-1)) {
        odlro.push_back(MC_Data::mean(n));
        odlro.push_back(MC_Data::stddev(n));
      }
      n += max_dist_;
    }

    /*
    for (int i=0; i<num_bond_types_; ++i) {
      for (int j=0; j<num_bond_types_; ++j) {
        if (d>=2) fs_ << MC_Data::result_str(n);
        // value at max distance
        if (d==(max_dist_-1)) {
          odlro.push_back(MC_Data::mean(n));
          odlro.push_back(MC_Data::stddev(n));
        }
        ++n;
      }
    }*/
    if (d>=2) {
      fs_ << MC_Data::conv_str(d); 
      fs_ << std::endl; 
    }
  }
  // Order parameter (from max_distance value)
  fs_ << "\n# phi(1)"; 
  for (const auto& p : xvals) fs_ << std::scientific << std::setw(14) << p;
  //fs_ << std::setw(6) << max_dist_-1; 
  for (int i=0; i<odlro.size(); i+=2) {
    fs_ << std::scientific << std::setw(15) << std::sqrt(std::abs(odlro[i]));
    fs_ << std::fixed << std::setw(10) << std::sqrt(std::abs(odlro[i+1]));
  }
  fs_ << MC_Data::conv_str(0); 
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


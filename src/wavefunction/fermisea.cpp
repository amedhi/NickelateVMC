/*---------------------------------------------------------------------------
* @Author: Amal Medhi, amedhi@mbpro
* @Date:   2019-02-20 12:21:42
* @Last Modified by:   amedhi
* @Last Modified time: 2020-04-03 11:29:40
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include <numeric>
#include "./fermisea.h"
#include <fstream>
#include <iomanip>

namespace var {

Fermisea::Fermisea(const input::Parameters& inputs, const lattice::LatticeGraph& graph) 
  : GroundState(true)
{
  init(inputs, graph);
}

int Fermisea::init(const input::Parameters& inputs, 
  const lattice::LatticeGraph& graph)
{
  name_ = "Fermisea";
  // sites & bonds
  num_sites_ = graph.num_sites();
  num_bonds_ = graph.num_bonds();
  // particle number
  set_particle_num(inputs);

  // build MF Hamiltonian
  varparms_.clear();
  mf_model_.init(graph.lattice());
  std::string name;
  double defval;
  using namespace model;
  model::CouplingConstant cc;

  if (graph.lattice().id()==lattice::lattice_id::SQUARE_NNN) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=1.0, inputs);
    cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"});
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
  }
  else if (graph.lattice().id()==lattice::lattice_id::SIMPLECUBIC) {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=1.0, inputs);
    mf_model_.add_parameter(name="th", defval=1.0, inputs);
    cc = CouplingConstant({0,"-t"},{1,"-t"},{2,"-tp"},{3,"-tp"},{4,"-th"});
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
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
    mf_model_.add_parameter(name="t_RN_200", defval=0.0, inputs);
    mf_model_.add_parameter(name="t_RN_202", defval=0.0, inputs);

    // bond operators
    cc.create(15);
    cc.add_type(0,  "t_NN_100");
    cc.add_type(1,  "t_NN_100");
    cc.add_type(2,  "t_NN_110");
    cc.add_type(3,  "t_NN_200");
    cc.add_type(4,  "t_NN_001");
    cc.add_type(5,  "t_RR_100");
    cc.add_type(6,  "t_RR_100");
    cc.add_type(7,  "t_RR_001");
    cc.add_type(8,  "t_RR_101");
    cc.add_type(9,  "t_RR_102");
    cc.add_type(10, "t_RR_110");
    cc.add_type(11, "t_RR_002");
    cc.add_type(12, "t_RR_111");
    cc.add_type(13, "t_RN_200");
    cc.add_type(14, "t_RN_202");
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());

    // site operators
    cc.create(2);
    cc.add_type(0, "e_N");
    cc.add_type(1, "e_R");
    mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());
  }
  else if (graph.lattice().id()==lattice::lattice_id::NICKELATE_2D) {
    mf_model_.add_parameter(name="e_N", defval=0.0, inputs);
    mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
    mf_model_.add_parameter(name="t_NN_100", defval=1.0, inputs);
    mf_model_.add_parameter(name="t_NN_110", defval=0.0, inputs);
    mf_model_.add_parameter(name="t_RR_100", defval=1.0, inputs);
    mf_model_.add_parameter(name="t_RR_110", defval=0.0, inputs);
    mf_model_.add_parameter(name="t_RN_200", defval=0.0, inputs);
    // bond operators
    cc.create(7);
    cc.add_type(0,  "t_NN_100");
    cc.add_type(1,  "t_NN_100");
    cc.add_type(2,  "t_NN_110");
    cc.add_type(3,  "t_RR_100");
    cc.add_type(4,  "t_RR_100");
    cc.add_type(5,  "t_RR_110");
    cc.add_type(6,  "t_RN_200");
    mf_model_.add_bondterm(name="hopping", cc, op::spin_hop());
    // site operators
    cc.create(2);
    cc.add_type(0, "e_N");
    cc.add_type(1, "e_R");
    mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());
  }
  else if (graph.lattice().id()==lattice::lattice_id::NICKELATE_2L) {
    mf_model_.add_parameter(name="e_R", defval=0.0, inputs);
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_parameter(name="tp", defval=0.0, inputs);
    mf_model_.add_parameter(name="th", defval=0.0, inputs);

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
    cc.add_type(0, "0");
    cc.add_type(1, "e_R");
    mf_model_.add_siteterm(name="ni_sigma", cc, op::ni_sigma());
  }
  else {
    mf_model_.add_parameter(name="t", defval=1.0, inputs);
    mf_model_.add_bondterm(name="hopping", cc="-t", op::spin_hop());
  }

  // chemical potential
  noninteracting_mu_ = true;
  // finalize MF Hamiltonian
  mf_model_.finalize(graph);
  num_varparms_ = varparms_.size();

  // bloch basis
  blochbasis_.construct(graph);
  num_kpoints_ = blochbasis_.num_kpoints();
  kblock_dim_ = blochbasis_.subspace_dimension();
  // FT matrix for transformation from 'site basis' to k-basis
  set_ft_matrix(graph);
  // work arrays
  work_.resize(kblock_dim_,kblock_dim_);
  phi_k_.resize(num_kpoints_);
  for (int k=0; k<num_kpoints_; ++k) {
    phi_k_[k].resize(kblock_dim_,kblock_dim_);
  } 

  // for calculating SC correlations 
  rmax_ =  graph.lattice().size1()/2+1;
  alpha_ = graph.lattice().basis_vector_a1();
  beta_ = graph.lattice().basis_vector_a2();
  R_list_.resize(rmax_);
  corr_aa_.resize(kblock_dim_,rmax_);
  corr_ab_.resize(kblock_dim_,rmax_);
  corr_fs_.resize(kblock_dim_,rmax_);
  for (int r=0; r<rmax_; ++r) {
    R_list_[r] = r*alpha_;
    //std::cout << r << "  " << R_list_[r].transpose() << "\n";
  }
  //getchar();

  return 0;
}

std::string Fermisea::info_str(void) const
{
  std::ostringstream info;
  info << "# Ground State: '"<<name_<<"'\n";
  info << "# Hole doping = "<<hole_doping()<<"\n";
  info << "# Particles = "<< num_upspins()+num_dnspins();
  info << " (Nup = "<<num_upspins()<<", Ndn="<<num_dnspins()<<")\n";
  info.precision(6);
  return info.str();
}

void Fermisea::update(const input::Parameters& inputs)
{
  // update from input parameters
  // hole doping might have changed
  set_particle_num(inputs);
  // update MF model
  mf_model_.update(inputs);
}

void Fermisea::update(const var::parm_vector& pvector, const unsigned& start_pos)
{
}

void Fermisea::get_wf_amplitudes(Matrix& psi) 
{
  construct_groundstate();
  get_pair_amplitudes(phi_k_);
  get_pair_amplitudes_sitebasis(phi_k_, psi);
  /*std::cout << "Fermisea::update: SC correlations\n";
  get_sc_correlation();
  getchar();
  */
}

void Fermisea::get_wf_gradient(std::vector<Matrix>& psi_gradient) 
{
  for (auto& mat : psi_gradient) mat.setZero();
}

void Fermisea::get_pair_amplitudes(std::vector<ComplexMatrix>& phi_k)
{
  for (int i=0; i<kshells_up_.size(); ++i) {
    int k = kshells_up_[i].k;
    int m = kshells_up_[i].nmax+1;
    Vector3d kvec = blochbasis_.kvector(k);
    mf_model_.construct_kspace_block(kvec);
    es_k_up.compute(mf_model_.quadratic_spinup_block());
    mf_model_.construct_kspace_block(-kvec);
    es_minusk_up.compute(mf_model_.quadratic_spinup_block());
    phi_k[k] = es_k_up.eigenvectors().block(0,0,kblock_dim_,m)
      		 * es_minusk_up.eigenvectors().transpose().block(0,0,m,kblock_dim_);
    /*
    std::cout << "kvec = "<< kvec.transpose() << "\n"; 
    std::cout << mf_model_.quadratic_spinup_block() << "\n"; 
    std::cout << phi_k[k] << "\n"; getchar();
    */
  }
}

void Fermisea::get_pair_amplitudes_sitebasis(const std::vector<ComplexMatrix>& phi_k, 
  Matrix& psi)
{
  // psi = FTU_ * PHI_K * conjugate(transpose(FTU_))
  // PHI_K is block diagonal (k-th block is phi_k) 
  int p = 0;
  for (int i=0; i<num_kpoints_; ++i) {
    int q = 0;
    for (int j=0; j<num_kpoints_; ++j) {
      work_.setZero();
      for (int ks=0; ks<kshells_up_.size(); ++ks) {
    	  int k = kshells_up_[ks].k;
        work_ += FTU_(i,k) * phi_k[k] * std::conj(FTU_(j,k));
      }
      // std::cout << work_ << "\n"; getchar();
      // copy transformed block
      //psi.block(p,q,kblock_dim_,kblock_dim_) = 
      for (int m=0; m<kblock_dim_; ++m) {
        for (int n=0; n<kblock_dim_; ++n) {
          psi(p+m,q+n) = ampl_part(work_(m,n));
      	}
      }
      q += kblock_dim_;
    }
    p += kblock_dim_;
  }
  //std::cout << psi << "\n";
  //getchar();
}

double Fermisea::get_mf_energy(void)
{
  double mf_energy = 0.0;
  return mf_energy/num_sites_;
}


void Fermisea::construct_groundstate(void)
{
  int num_upspins_ = num_upspins();
  int num_dnspins_ = num_dnspins();
  if (have_TP_symmetry_ || num_upspins_!=num_dnspins_) {
    /* Has T.P (Time Reversal * Inversion) symmetry. 
       So we have e_k(up) = e_k(dn).
    */
    std::vector<std::pair<int,int>> qn_list; // list of (k,n)
    std::vector<double> ek;
    for (int k=0; k<num_kpoints_; ++k) {
      Vector3d kvec = blochbasis_.kvector(k);
      mf_model_.construct_kspace_block(kvec);
      es_k_up.compute(mf_model_.quadratic_spinup_block(), Eigen::EigenvaluesOnly);
      ek.insert(ek.end(),es_k_up.eigenvalues().data(),
        es_k_up.eigenvalues().data()+kblock_dim_);
      //std::cout << kvec.transpose() << " " << es_k_up.eigenvalues() << "\n"; getchar();
      for (int n=0; n<kblock_dim_; ++n) {
        qn_list.push_back({k, n});
      }
    }
    //std::sort(ek.begin(),ek.end());
    //for (const auto& e : ek) std::cout << e << "\n";

    // Indices in the original ek-array is sorted according to increasing ek
    std::vector<int> idx(ek.size());
    std::iota(idx.begin(),idx.end(),0);
    std::sort(idx.begin(),idx.end(),[&ek](const int& i1, const int& i2) 
      {return ek[i1]<ek[i2];});
    /*for (int i=0; i<ek.size(); ++i) {
      //std::cout << i << "  " << idx[i] << "  " << ek[idx[i]] << "\n";
      std::cout << i+1 << " " << blochbasis_.kvector(idx[i]).transpose()  << "\n";
    }
    getchar();
    */
    
    // mean energy
    /*double e0 = 0.0;
    for (int i=0; i<num_upspins_; ++i) {
      e0 += ek[idx[i]];
    }
    e0 = 2 * e0 / num_sites_;
    std::cout << "e0 = " << e0 << "\n"; getchar();
    */

    // check for degeneracy 
    double degeneracy_tol = 1.0E-12;
    int top_filled_level = num_upspins()-1;
    int num_degen_states = 1;
    int num_valence_particle = 1;
    fermi_energy_ = ek[idx[top_filled_level]];
    total_energy_ = 0.0;
    for (int i=0; i<=top_filled_level; ++i) {
      total_energy_ += ek[idx[i]];
    }
    total_energy_ = 2.0*total_energy_/num_sites_;
    std::cout << "Total KE = " << total_energy_ << "\n";

    // look upward in energy
    for (int i=top_filled_level+1; i<ek.size(); ++i) {
      if (std::abs(fermi_energy_-ek[idx[i]])>degeneracy_tol) break;
      num_degen_states++;
    }
    // look downward in energy
    if (num_degen_states>1) {
      for (int i=top_filled_level-1; i>=0; --i) {
        if (std::abs(fermi_energy_ - ek[idx[i]])>degeneracy_tol) break;
        num_degen_states++;
        num_valence_particle++;
      }
      // warn
      //if (degeneracy_warning_) {
        std::cout << ">> warning: groundstate degeneracy: " << 2*num_valence_particle <<
          " particles in " << 2*num_degen_states << " states." << "\n";
      //}
    }
    /* 
      Filled k-shells. A k-shell is a group of energy levels having same 
      value of quantum number k.
    */
    // find 'nmax' values of filled k-shells
    std::vector<int> shell_nmax(num_kpoints_);
    for (auto& elem : shell_nmax) elem = -1; // invalid default value
    int k, n;
    for (int i=0; i<num_upspins_; ++i) {
      int state = idx[i]; 
      std::tie(k,n) = qn_list[state];
      if (shell_nmax[k] < n) shell_nmax[k] = n;
    }
    // store the filled k-shells
    kshells_up_.clear();
    for (int k=0; k<num_kpoints_; ++k) {
      int nmax = shell_nmax[k];
      if (nmax != -1) kshells_up_.push_back({k,0,nmax});
    }
    
    /*
    for (unsigned k=0; k<kshells_up_.size(); ++k) {
      std::cout << kshells_up_[k].k << " " << kshells_up_[k].nmin << "  "
          << kshells_up_[k].nmax << "\n";
    }
    getchar();
    */
    
  }
  else {
    /* 
      No T.P (Time Reversal * Inversion) symmetry. 
      Need to consider both UP and DN states.
    */
    throw std::range_error("Fermisea::construct_groundstate: case not implemented\n");
  }
}

void Fermisea::get_sc_correlation(void)
{
  corr_aa_.setZero();
  corr_ab_.setZero();
  corr_fs_.setZero();
  for (int r=1; r<rmax_; ++r) {
    Vector3d R = R_list_[r];
    for (int i=0; i<kshells_up_.size(); ++i) {
      int k = kshells_up_[i].k;
      Vector3d kvec = blochbasis_.kvector(k);
      auto k_ab = std::exp(-ii()*kvec.dot(alpha_-beta_));
      // 1-body correlations
      for (int n=0; n<kblock_dim_; ++n) {
        corr_fs_(n,r) += std::exp(-ii()*kvec.dot(R));
      }
      for (int j=0; j<kshells_up_.size(); ++j) {
        int kp = kshells_up_[j].k;
        Vector3d kpvec = blochbasis_.kvector(kp);
        Vector3d kkplus = kvec+kpvec;
        Vector3d kkminus = kvec-kpvec;
        auto kkp = std::exp(ii()*kkplus.dot(R));
        auto kp_ab = std::exp(ii()*kpvec.dot(beta_-alpha_));
        auto term1 = 2*std::cos(kkminus.dot(alpha_));
        auto term2 = std::exp(ii()*(-kvec.dot(alpha_)+kpvec.dot(beta_)))
                    +std::exp(ii()*(kvec.dot(beta_)-kpvec.dot(alpha_)));
        for (int n=0; n<kblock_dim_; ++n) {
          if (n<=kshells_up_[i].nmax && n<=kshells_up_[j].nmax) {
            corr_aa_(n,r) += kkp*(2.0 + term1);
            corr_ab_(n,r) += kkp*(k_ab + term2 + kp_ab);
          }
        }
      }
    }
    for (int n=0; n<kblock_dim_; ++n) {
      corr_aa_(n,r) /= (2*num_kpoints_*num_kpoints_);
      corr_ab_(n,r) /= (2*num_kpoints_*num_kpoints_);
      corr_fs_(n,r) /= num_kpoints_;
    }
  }

  // print 
  std::ofstream fs("fs_sccorr.txt");
  fs << std::right<<std::scientific<<std::uppercase<<std::setprecision(6);
  for (int r=1; r<rmax_; ++r) {
    fs<<std::setw(6)<<r<< " ";
    for (int n=0; n<kblock_dim_; ++n) {
      //fs<<std::setw(15)<<std::real(corr_fs_(n,r)); 
      //fs<<std::setw(15)<<std::imag(corr_fs_(n,r)); 
      fs<<std::setw(15)<<std::real(corr_aa_(n,r)); 
      fs<<std::setw(15)<<std::imag(corr_aa_(n,r)); 
      fs<<std::setw(15)<<std::real(corr_ab_(n,r)); 
      fs<<std::setw(15)<<std::imag(corr_ab_(n,r)); 
    }
    fs<<"\n";
  }
  fs.close();
  std::cout << "sc_correlation done\n";
  //getchar();
}


} // end namespace var

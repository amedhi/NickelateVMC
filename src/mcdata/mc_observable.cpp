/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-06 11:31:00
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-16 09:32:43
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#include "./mc_observable.h"
#include <boost/algorithm/string.hpp>

namespace mcdata {

MC_Observable::MC_Observable() 
  : MC_Data() 
{
  num_dataset_=0;
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
}

MC_Observable::MC_Observable(const std::string& name, const int& size,
  const bool& no_error_bar, const bool& replace_mode) 
{
  this->init(name,size,no_error_bar);
  replace_mode_ = replace_mode;
} 

void MC_Observable::init(const std::string& name, const int& size, const bool& no_error_bar)
{
  MC_Data::init(name,size,no_error_bar);
  name_ = name;
  num_dataset_=0;
  avg_mcdata_.init(name,size,no_error_bar);
  avg_stddev_.resize(size);
  avg_stddev_.setZero();
  avg_tau_ = 0.0;
  //elem_names_ = {name};
  elem_names_.resize(size);
  // default file name
  fname_ = name_;
  boost::to_lower(fname_);
  auto pos = fname_.find('^');
  if (pos != std::string::npos) fname_.erase(pos,1);
  fname_ = "res_"+fname_+".txt";
}

void MC_Observable::set_ofstream(const std::string& prefix)
{
  // file name
  fname_ = name_;
  boost::to_lower(fname_);
  auto pos = fname_.find('^');
  if (pos != std::string::npos) fname_.erase(pos,1);
  fname_ = prefix+"res_"+fname_+".txt";
}

void MC_Observable::resize(const int& size)
{
  MC_Data::resize(size);
  avg_mcdata_.resize(size);
  avg_stddev_.resize(size);
  avg_stddev_.setZero(size);
  elem_names_.resize(size);
  for (unsigned i=0; i<size; ++i) 
    elem_names_[i] = "elem" + std::to_string(i);
}

void MC_Observable::resize(const int& size, const std::vector<std::string>& elem_names)
{
  MC_Data::resize(size);
  avg_mcdata_.resize(size);
  avg_stddev_.resize(size);
  avg_stddev_.setZero(size);
  elem_names_ = elem_names;
}

void MC_Observable::check_on(const input::Parameters& inputs, const bool& replace_mode) 
{
  int no_warn;
  is_on_ = inputs.set_value(name(), false, no_warn);
  replace_mode_ = replace_mode;
}

void MC_Observable::print_heading(const std::string& header,
  const std::vector<std::string>& xvars) 
{
  if (!is_on()) return;
  if (heading_printed_) return;
  if (!replace_mode_) return;
  if (!fs_.is_open()) open_file();
  fs_ << header;
  fs_ << "# Results: " << name() << "\n";
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# ";
  fs_ << std::left;
  //fs_ << std::setw(14)<<xvar_name;
  for (const auto& p : xvars) fs_ << std::setw(15)<<p.substr(0,15);
  // total value
  if (MC_Data::size()>1 && have_total_)
    fs_ << std::setw(14)<<"Total"<<std::setw(11)<<"err";
  for (const auto& name : elem_names_) 
    fs_ << std::setw(14)<<name<<std::setw(11)<<"err";
  //fs_ << std::setw(9)<<"samples";
  fs_ << std::setw(9)<<"samples"<<std::setw(12)<<"converged"<<std::setw(6)<<"tau";
  fs_ << std::endl;
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << std::flush;
  heading_printed_ = true;
  close_file();
}

void MC_Observable::print_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return;
  if (!fs_.is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  // total value
  //if (MC_Data::size()>1)
  if (MC_Data::size()>1 && have_total_)
    fs_ << MC_Data::result_str(-1); 
  for (int i=0; i<MC_Data::size(); ++i) 
    fs_ << MC_Data::result_str(i); 
  fs_ << MC_Data::conv_str(0); //.substr(0,10); 
  fs_ << std::endl;
  fs_ << std::flush;
  close_file();
} 

void MC_Observable::print_grand_result(const std::vector<double>& xpvals) 
{
  if (!is_on()) return;
  if (!fs_.is_open()) open_file();
  fs_ << std::right;
  fs_ << std::scientific << std::uppercase << std::setprecision(6);
  fs_ << "#" << std::string(72, '-') << "\n";
  fs_ << "# grand average:\n"; // 
  fs_ << "#" << std::string(72, '-') << "\n";
  for (const auto& p : xpvals) 
    fs_ << std::setw(14) << p;
  // total value
  if (avg_mcdata_.size()>1 && have_total_)
    fs_ << avg_mcdata_.result_str(-1); 
  for (int i=0; i<avg_mcdata_.size(); ++i) {
    fs_ << std::scientific << std::uppercase << std::setprecision(6) << std::right;
    fs_ << std::setw(15) << avg_mcdata_.mean(i); 
    fs_ << std::fixed << std::setw(10) << avg_stddev_[i]/num_dataset_; 
  }
  fs_ << std::endl;
  fs_ << std::flush;
  close_file();
} 

void MC_Observable::open_file(void) 
{
  if (fs_.is_open()) return;
  if (replace_mode_) {
    fs_.open(fname_);
    replace_mode_ = false;
  }
  else fs_.open(fname_, std::ios::app);
  if (!fs_.is_open()) 
    throw std::runtime_error("Observable::open_file: file open failed");
}

void MC_Observable::close_file(void) 
{
  fs_.close();
}

void MC_Observable::reset_grand_data(void) 
{ 
  num_dataset_=0; 
  avg_mcdata_.clear(); 
  avg_stddev_.setZero(); 
  avg_tau_ = -1;
} 

void MC_Observable::save_result(void)
{
  avg_mcdata_ << MC_Data::mean_data();
  avg_stddev_ += MC_Data::stddev_data();
  //std::cout << "mean_data="<< MC_Data::mean_data() << "\n"; 
  //std::cout << "save_result="<< avg_mcdata_.result_str(-1) << "\n"; 
  avg_tau_ += MC_Data::tau();
  ++num_dataset_;
}

void MC_Observable::avg_grand_data(void) 
{ 
  //std::cout << avg_mcdata_.result_str(-1); 
  mcdata::MC_Data::copy_finalize(avg_mcdata_);
  //std::cout << MC_Data::result_str(-1); 
  //getchar();
} 

MC_Data& MC_Observable::grand_data(void) 
{
  return avg_mcdata_; 
} 

data_t MC_Observable::grand_stddev(void) const 
{ 
  if (num_dataset_ > 0) return avg_stddev_/static_cast<double>(num_dataset_); 
  else return avg_stddev_;
}


} // end namespace mcdata



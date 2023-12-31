/*---------------------------------------------------------------------------
* Author: Amal Medhi
* Date:   2017-05-06 11:14:26
* Last Modified by:   Amal Medhi, amedhi@macbook
* Last Modified time: 2017-05-16 09:32:41
* Copyright (C) Amal Medhi, amedhi@iisertvm.ac.in
*----------------------------------------------------------------------------*/
#ifndef MC_OBSERVABLE_H
#define MC_OBSERVABLE_H

#include <string>
#include <vector>
#include <map>
#include <stdexcept>
#include <fstream>
#include <iomanip>
#include "../scheduler/task.h"
#include "../mcdata/mcdata.h"

namespace mcdata {

class MC_Observable : public mcdata::MC_Data
{
public:
  using data_t = mcdata::data_t;
  using scalar_t = mcdata::scalardata_t;
  MC_Observable();
  MC_Observable(const std::string& name, const int& size=1, const bool& no_error_bar=false,
    const bool& replace_mode=true);
  ~MC_Observable() {}
  void init(const std::string& name, const int& size=1, const bool& no_error_bar=false) override; 
  void set_ofstream(const std::string& prefix);
  void resize(const int& size) override;
  void resize(const int& size, const std::vector<std::string>& elem_names);
  void set_replace_mode(const bool& replace_mode) { replace_mode_=replace_mode; }
  void set_have_total(void) { have_total_=true; }
  virtual void reset(void) { MC_Data::clear(); }
  void check_on(const input::Parameters& inputs, const bool& replace_mode); 
  void switch_on(void) { is_on_=true; }
  void switch_off(void) { is_on_=false; if (fs_.is_open()) fs_.close(); }
  operator int(void) const { return is_on(); }
  const bool& is_on(void) const { return is_on_; }
  void open_file(void); 
  void close_file(void); 
  bool is_open(void) const { return fs_.is_open(); }
  std::ofstream& fs(void) { return fs_; }  
  MC_Data& grand_data(void);
  data_t grand_stddev(void) const;
  void reset_grand_data(void); 
  void save_result(void);
  void avg_grand_data(void); 
  virtual void print_heading(const std::string& header, 
    const std::vector<std::string>& xvars);
  virtual void print_result(const std::vector<double>& xvals); 
  virtual void print_grand_result(const std::vector<double>& xvals); 
protected:
  // for printing
  std::vector<std::string> elem_names_;
  std::ofstream fs_;
  bool heading_printed_{false};
  bool replace_mode_{true};
private:
  int num_dataset_{0};
  MC_Data avg_mcdata_;
  data_t avg_stddev_;
  double avg_tau_;
  bool is_on_{false};
  bool have_total_{false};
  std::string name_{""};
  std::string fname_{""};
  mcdata::MC_Data data(void) { return mcdata::MC_Data(*this); }
};


} // end namespace mcdata

#endif
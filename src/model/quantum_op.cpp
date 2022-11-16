/*---------------------------------------------------------------------------
* Copyright (C) 2015-2016 by Amal Medhi <amedhi@iisertvm.ac.in>.
* @Author: Amal Medhi
* @Date:   2022-11-13 21:53:39
* @Last Modified by:   Amal Medhi
* @Last Modified time: 2022-11-14 21:46:41
*----------------------------------------------------------------------------*/
#include <set>
#include "./quantum_op.h"

namespace model {

ProjectionOp::ProjectionOp(const projection_t& pj)
{
  this->clear();
  // expr is applicable for all site types
  super_type::insert({0, pj});
  global_type_ = true;
}

ProjectionOp::ProjectionOp(const value_type& type0, const value_type& type1, 
	const value_type& type2, const value_type& type3, const value_type& type4,
  const value_type& type5)
{
  set(type0, type1, type2, type3, type4, type5);
}

void ProjectionOp::clear(void) 
{
	super_type::clear();	
  is_present_ = false;
  global_type_ = false;
  finalized_ = false;
}

ProjectionOp& ProjectionOp::operator=(const projection_t& pj)
{
  this->clear();
  // expr is applicable for all site types
  super_type::insert({0, pj});
  global_type_ = true;
  return *this;
}

void ProjectionOp::set(const value_type& type0, const value_type& type1, 
	const value_type& type2, const value_type& type3, const value_type& type4, 
  const value_type& type5)
{
  this->clear();
  super_type::insert(type0);
  if (type1.second != pjn::NONE) {
    super_type::insert(type1); 
  }
  if (type2.second != pjn::NONE) {
    super_type::insert(type2); 
  }
  if (type3.second != pjn::NONE) {
    super_type::insert(type3); 
  }
  if (type4.second != pjn::NONE) {
    super_type::insert(type4); 
  }
  if (type5.second != pjn::NONE) {
    super_type::insert(type5); 
  }
  global_type_ = false;
}

void ProjectionOp::finalize(const unsigned& num_site_types)
{
	if (super_type::size() == 0) {
  	super_type::insert({0,projection_t::NONE});
  	is_present_ = false;
  	global_type_ = true;
	}

	if (global_type_) {
		if (super_type::at(0)==projection_t::NONE) {
  		is_present_ = false;
		}
		else {
  		is_present_ = true;
		}
	}
	else {
		// check if size match
		if (super_type::size() != num_site_types) {
    	throw std::range_error("model::ProjectionOp::finalize: size mismatch");
		}
		// check if all types are defined
		is_present_ = false;
		for (unsigned i=0; i<num_site_types; ++i) {
			if (super_type::find(i) == super_type::end()) {
    		throw std::range_error("model::ProjectionOp::finalize: invalid definition");
			}
			if (super_type::at(i)!=projection_t::NONE) {
				is_present_ = true;
			}
		} 
	}
	finalized_ = true;
}

const projection_t& ProjectionOp::get(const unsigned& i) const
{
	if (!finalized_) {
    throw std::logic_error("model::ProjectionOp::get: operator not finalized");
	}

	if (global_type_) {
		return super_type::at(0);
	}

	if (super_type::find(i)==super_type::end()) {
    throw std::range_error("model::ProjectionOp::get: operator not found");
	}
	return super_type::at(i);
}


} // namespace model
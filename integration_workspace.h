// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2011-2017 Free Software Foundation, Inc.
//
// This file is part of the GNU ISO C++ Library.  This library is free
// software; you can redistribute it and/or modify it under the
// terms of the GNU General Public License as published by the
// Free Software Foundation; either version 3, or (at your option)
// any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along
// with this library; see the file COPYING3.  If not see
// <http://www.gnu.org/licenses/>.
//
// Implements the integration_workspace class which stores temporary data
// for performing integrals
// Based on gsl/integration/workspace.c

#ifndef INTEGRATION_WORKSPACE_H
#define INTEGRATION_WORKSPACE_H 1

#include <vector>
#include <limits>
#include <cmath>

namespace __gnu_test
{

  template<typename _Tp>
    class integration_workspace
    {
    private:

      struct interval
      {
	_Tp _M_lower_lim;
	_Tp _M_upper_lim;
	_Tp _M_result;
	_Tp _M_abs_error;
	std::size_t _M_order;
	std::size_t _M_level;
      };

      /**
       * Comparison of cquad intervals.
       */
      struct interval_comp
      {
	bool
	operator()(const interval& __ivl,
		   const interval& __ivr)
	{ return __ivl._M_abs_error < __ivr._M_abs_error; }
      };

      std::size_t _M_capacity;
      std::size_t _M_size;
      std::size_t _M_maxerr_index;
      std::size_t _M_maximum_level;
      std::vector<interval> _M_ival;
      bool _M_try_resize;

    public:

      integration_workspace(std::size_t __cap)
      : _M_capacity(__cap),
	_M_size(0),
	_M_maxerr_index(0),
	_M_maximum_level(0),
	_M_ival(__cap)
      { }

      void
      resize()
      {
	const auto __new_cap = std::size_t(1 + 1.5 * this->capacity());
	this->_M_capacity = __new_cap;
	this->_M_ival.resize(__new_cap);
      }

      void
      set_initial_limits(_Tp a0, _Tp b0)
      {
	this->_M_size = 0;
	this->_M_maxerr_index = 0;
	if (this->_M_capacity > 0)
	  {
	    this->_M_ival[0]._M_lower_lim = a0;
	    this->_M_ival[0]._M_upper_lim = b0;
	    this->_M_ival[0]._M_result = _Tp{0};
	    this->_M_ival[0]._M_abs_error = _Tp{0};
	    this->_M_ival[0]._M_order = 0;
	    this->_M_ival[0]._M_level = 0;
	  }
	this->_M_maximum_level = 0;
      }

      void
      set_initial_results(_Tp result, _Tp error)
      {
	this->_M_size = 1;
	this->_M_maxerr_index = 0;
	if (this->_M_capacity > 0)
	  {
	    this->_M_ival[0]._M_result = result;
	    this->_M_ival[0]._M_abs_error = error;
	    this->_M_ival[0]._M_order = 0;
	    this->_M_ival[0]._M_level = 0;
	  }
	this->_M_maximum_level = 0;
      }

      void
      set_initial(_Tp __a0, _Tp __b0, _Tp __result0, _Tp __error0)
      {
	this->_M_size = 1;
	this->_M_maxerr_index = 0;
	if (this->_M_capacity > 0)
	  {
	    this->_M_ival[0]._M_lower_lim = __a0;
	    this->_M_ival[0]._M_upper_lim = __b0;
	    this->_M_ival[0]._M_result = __result0;
	    this->_M_ival[0]._M_abs_error = __error0;
	    this->_M_ival[0]._M_order = 0;
	    this->_M_ival[0]._M_level = 0;
	  }
	this->_M_maximum_level = 0;
      }

      void sort_error();
      void sort_results();

      void append(_Tp __a, _Tp __b, _Tp __area, _Tp __error);

      void update(_Tp __a1, _Tp __b1, _Tp __area1, _Tp __error1,
		  _Tp __a2, _Tp __b2, _Tp __area2, _Tp __error2);

      void
      retrieve(_Tp& __lolim, _Tp& __uplim, _Tp& __res, _Tp& __err) const
      {
	const auto __ii = this->maxerr_order();
	__lolim = this->_M_ival[__ii]._M_lower_lim;
	__uplim = this->_M_ival[__ii]._M_upper_lim;
	__res = this->_M_ival[__ii]._M_result;
	__err = this->_M_ival[__ii]._M_abs_error;
      }

      std::size_t
      size() const
      { return this->_M_size; }

      std::size_t
      capacity() const
      { return this->_M_capacity; }

      std::size_t
      maxerr_order() const
      { return this->_M_ival[this->_M_maxerr_index]._M_order; }

      _Tp
      lower_lim(std::size_t __ii) const
      { return this->_M_ival[__ii]._M_lower_lim; }

      _Tp
      upper_lim(std::size_t __ii) const
      { return this->_M_ival[__ii]._M_upper_lim; }

      _Tp
      result(std::size_t __ii) const
      { return this->_M_ival[__ii]._M_result; }

      _Tp
      abs_error(std::size_t __ii) const
      { return this->_M_ival[__ii]._M_abs_error; }

      size_t
      order(std::size_t __ii) const
      { return this->_M_ival[__ii]._M_order; }

      size_t
      level(std::size_t ii) const
      { return this->_M_ival[ii]._M_level; }

      // Only used by qagp
      _Tp
      set_abs_error(std::size_t __ii, _Tp __abserr)
      { return this->_M_ival[__ii]._M_abs_error = __abserr; }

      // Only used by qagp
      void
      set_level(std::size_t __ii, std::size_t __lvl)
      { this->_M_ival[__ii]._M_level = __lvl; }

      std::size_t
      current_level() const
      { return this->_M_ival[this->maxerr_order()]._M_level; }

      std::size_t
      max_level() const
      { return this->_M_maximum_level; }

      bool
      large_interval() const
      {
	if (this->_M_ival[this->maxerr_order()]._M_level
		 < this->_M_maximum_level)
	  return true;
	else
	  return false;
      }

      bool
      increase_maxerr_index()
      {
	if (this->_M_size < 2)
	  return false;

	auto __id = this->_M_maxerr_index;
	auto __last = this->_M_size - 1;

	std::size_t __jupbnd;
	if (__last > (1 + this->_M_capacity / 2))
	  __jupbnd = this->_M_capacity + 1 - __last;
	else
	  __jupbnd = __last;

	for (auto __k = __id; __k <= __jupbnd; ++__k)
	  {
	    auto __i_max = this->_M_ival[this->_M_maxerr_index]._M_order;
	    if (this->_M_ival[__i_max]._M_level < this->_M_maximum_level)
	      return true;
	    ++this->_M_maxerr_index;
	  }
	return false;
      }

      void
      reset_maxerr_index()
      { this->_M_maxerr_index = 0; }

      std::size_t
      get_maxerr_index() const
      { return this->_M_maxerr_index; }

      void
      set_maxerr_index(std::size_t maxerr_index)
      { this->_M_maxerr_index = maxerr_index; }

      _Tp
      sum_results() const
      {
	auto __result_sum = _Tp{0};
	for (std::size_t __kk = 0; __kk < this->_M_size; ++__kk)
	  __result_sum += this->_M_ival[__kk]._M_result;

	return __result_sum;
      }

      static bool
      subinterval_too_small(_Tp __a1, _Tp __a2, _Tp __b2)
      {
	const auto _S_eps = 100 * std::numeric_limits<_Tp>::epsilon();
	const auto _S_min = 1000 * std::numeric_limits<_Tp>::min();

	_Tp __tmp = (1 + _S_eps) * (std::abs(__a2) + _S_min);

	return std::abs(__a1) <= __tmp
	    && std::abs(__b2) <= __tmp;
      }
    };

} // namespace __gnu_test

#include "integration_workspace.tcc"

#endif // INTEGRATION_WORKSPACE_H

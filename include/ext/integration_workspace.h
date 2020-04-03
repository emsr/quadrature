// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
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

#include <algorithm>
#include <vector>
#include <limits>
#include <cmath>
#include <iosfwd>

namespace __gnu_cxx
{

  template<typename _Tp, typename _RetTp>
    class integration_workspace
    {
    private:

      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _ErrorTp = decltype(std::abs(_AreaTp{}));

      struct interval
      {
	_Tp __lower_lim;
	_Tp __upper_lim;
	_AreaTp __result;
	_ErrorTp __abs_error;
	std::size_t __depth;

	bool
	operator<(const interval& __iv)
	{ return this->__abs_error < __iv.__abs_error; }
      };

      /**
       * Comparison of quadrature intervals.
       * The comparison is by absolute error.
       */
      struct interval_comp
      {
	bool
	operator()(const interval& __ivl,
		   const interval& __ivr)
	{ return __ivl.__abs_error < __ivr.__abs_error; }
      };

      // The start of the heap.
      // This allows to skip the actual max error.
      std::size_t _M_curr_index;

      // The current maximum depth.
      std::size_t _M_max_depth;

      // The maximum size of the workspace.
      std::size_t _M_max_size;

      std::vector<interval> _M_ival;

    public:

      integration_workspace(std::size_t __cap)
      : _M_curr_index{0},
	_M_max_depth{0},
	_M_max_size(__cap),
	_M_ival{}
      {
	_M_ival.reserve(__cap);
      }

      void sort_error();

      void append(_Tp __a, _Tp __b, _AreaTp __area, _ErrorTp __error,
		  std::size_t __depth = 0);

      void split(_Tp __ab,
		 _AreaTp __area1, _ErrorTp __error1,
		 _AreaTp __area2, _ErrorTp __error2);

      const interval&
      retrieve() const
      { return this->_M_ival[this->curr_index()]; }

      std::size_t
      size() const
      { return this->_M_ival.size(); }

      std::size_t
      max_size() const
      { return this->_M_max_size; }

      std::size_t
      capacity() const
      { return this->_M_ival.capacity(); }

      /**
       * Return the current segment - the beginning of the vector + start
       * as an iterator.
       */
      typename std::vector<interval>::iterator
      begin()
      { return this->_M_ival.begin() + this->curr_index(); }

      /**
       * Return the end iterator of the vector.
       */
      typename std::vector<interval>::iterator
      end()
      { return this->_M_ival.end(); }

      /**
       * Return the current segment - the top of the heap - as an lvalue.
       */
      const interval&
      top() const
      { return this->_M_ival[this->curr_index()]; }

      /**
       * Return the current segment - the top of the heap - as an rvalue.
       */
      interval&
      top()
      { return this->_M_ival[this->curr_index()]; }

      void
      clear()
      {
	this->_M_curr_index = 0;
	this->_M_max_depth = 0;
	this->_M_ival.clear();
      }

      /**
       * Push a new segment into the current heap.
       * N.B. this->begin() includes curr_index()!
       */
      void
      push(const interval& __iv)
      {
	this->_M_ival.push_back(__iv);
	std::push_heap(this->begin(), this->end(), interval_comp{});
      }

      /**
       * Pop the current segment out of teh heap and erase it from the vector.
       * N.B. this->begin() includes curr_index()!
       */
      void
      pop()
      {
	std::pop_heap(this->begin(), this->end());
	this->_M_ival.pop_back();
      }

      /**
       * Return the lower limit for the segment at start + ii.
       */
      _Tp
      lower_lim(std::size_t __ii = 0) const
      { return this->_M_ival[this->curr_index() + __ii].__lower_lim; }

      /**
       * Return the upper limit for the segment at start + ii.
       */
      _Tp
      upper_lim(std::size_t __ii = 0) const
      { return this->_M_ival[this->curr_index() + __ii].__upper_lim; }

      /**
       * Return the integration result for the segment at start + ii.
       */
      _AreaTp
      result(std::size_t __ii = 0) const
      { return this->_M_ival[this->curr_index() + __ii].__result; }

      /**
       * Return the absolute error for the segment at start + ii.
       */
      _ErrorTp
      abs_error(std::size_t __ii = 0) const
      { return this->_M_ival[this->curr_index() + __ii].__abs_error; }

      /**
       * Return the subdivision depth for the segment at start + ii.
       */
      size_t
      depth(std::size_t __ii = 0) const
      { return this->_M_ival[this->curr_index() + __ii].__depth; }

      /**
       * Set the absolute error of the segment at start + ii to the given value.
       * Only used by qagp.
       */
      _ErrorTp
      set_abs_error(std::size_t __ii, _ErrorTp __abserr)
      { return this->_M_ival[this->curr_index() + __ii].__abs_error = __abserr; }

      /**
       * Set the depth of segment at start + ii to the given value.
       * Only used by qagp.
       */
      void
      set_depth(std::size_t __ii, std::size_t __d)
      { this->_M_ival[this->curr_index() + __ii].__depth = __d; }

      bool increment_curr_index();

      /**
       * Reset the index of the current segment to zero.
       * The entire segment vector will me made into a heap.
       */
      void
      reset_curr_index()
      {
	this->_M_curr_index = 0;
	this->sort_error();
      }

      /**
       * Return the depth of the current segment.
       */
      std::size_t
      curr_depth() const
      { return this->_M_ival[this->curr_index()].__depth; }

      std::size_t
      max_depth() const
      { return this->_M_max_depth; }

      /**
       * Return the index of the current segment in the segment vector.
       * The segments in [start, end) will be maintained in a heap.
       */
      std::size_t
      curr_index() const
      { return this->_M_curr_index; }

      /**
       * Return true if the subdivision depth of the current segment
       * is less than the maximum depth.
       */
      bool
      large_interval() const
      {
	if (this->curr_depth() < this->max_depth())
	  return true;
	else
	  return false;
      }

      /// Return the total integral:
      /// the sum of the results over all integration segments.
      _AreaTp
      total_integral() const
      {
	auto __result_sum = _AreaTp{0};
	for (const auto& __iv : this->_M_ival)
	  __result_sum += __iv.__result;
	return __result_sum;
      }

      /// Return the sum of the absolute errors over all integration segments.
      _ErrorTp
      total_error() const
      {
	auto __tot_error = _ErrorTp{0};
	for (auto& __iv : this->_M_ival)
	  __tot_error += __iv.__abs_error;
	return __tot_error;
      }

      /// Return the vector of integration intervals.
      const std::vector<interval>&
      intervals() const
      { return this->_M_ival; }

      /*
       * 
       */
      static bool
      subinterval_too_small(_Tp __a1, _Tp __a2, _Tp __b2)
      {
	const auto _S_eps = _Tp{100} * std::numeric_limits<_Tp>::epsilon();
	const auto _S_min = _Tp{1000} * std::numeric_limits<_Tp>::min();

	const auto __tmp = (_Tp{1} + _S_eps) * (std::abs(__a2) + _S_min);

	return std::abs(__a1) <= __tmp
	    && std::abs(__b2) <= __tmp;
      }
    };

  template<typename _Tp, typename _RetTp>
    std::ostream&
    operator<<(std::ostream& __out,
	       const integration_workspace<_Tp, _RetTp>& __ws);

} // namespace __gnu_cxx

#include <ext/integration_workspace.tcc>

#endif // INTEGRATION_WORKSPACE_H

//
// Copyright (C) 2011-2020 Free Software Foundation, Inc.
// Copyright (C) 2021-2022 Edward M. Smith-Rowland
//
// This program is free software; you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation; either version 3 of the License, or (at
// your option) any later version.
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

namespace emsr
{

  template<typename Tp, typename RetTp>
    class integration_workspace
    {
    private:

      using AreaTp = decltype(RetTp{} * Tp{});
      using ErrorTp = decltype(std::abs(AreaTp{}));

      struct interval
      {
	Tp lower_lim;
	Tp upper_lim;
	AreaTp result;
	ErrorTp abs_error;
	std::size_t depth;

	bool
	operator<(const interval& iv)
	{ return this->abs_error < iv.abs_error; }
      };

      /**
       * Comparison of quadrature intervals.
       * The comparison is by absolute error.
       */
      struct interval_comp
      {
	bool
	operator()(const interval& ivl,
		   const interval& ivr)
	{ return ivl.abs_error < ivr.abs_error; }
      };

      // The start of the heap.
      // This allows to skip the actual max error.
      std::size_t m_curr_index;

      // The current maximum depth.
      std::size_t m_max_depth;

      // The maximum size of the workspace.
      std::size_t m_max_size;

      std::vector<interval> m_ival;

    public:

      integration_workspace(std::size_t cap)
      : m_curr_index{0},
	m_max_depth{0},
	m_max_size(cap),
	m_ival{}
      {
	m_ival.reserve(cap);
      }

      void sort_error();

      void append(Tp a, Tp b, AreaTp area, ErrorTp error,
		  std::size_t depth = 0);

      void split(Tp ab,
		 AreaTp area1, ErrorTp error1,
		 AreaTp area2, ErrorTp error2);

      const interval&
      retrieve() const
      { return this->m_ival[this->curr_index()]; }

      std::size_t
      size() const
      { return this->m_ival.size(); }

      std::size_t
      max_size() const
      { return this->m_max_size; }

      std::size_t
      capacity() const
      { return this->m_ival.capacity(); }

      /**
       * Return the current segment - the beginning of the vector + start
       * as an iterator.
       */
      typename std::vector<interval>::iterator
      begin()
      { return this->m_ival.begin() + this->curr_index(); }

      /**
       * Return the end iterator of the vector.
       */
      typename std::vector<interval>::iterator
      end()
      { return this->m_ival.end(); }

      /**
       * Return the current segment - the top of the heap - as an lvalue.
       */
      const interval&
      top() const
      { return this->m_ival[this->curr_index()]; }

      /**
       * Return the current segment - the top of the heap - as an rvalue.
       */
      interval&
      top()
      { return this->m_ival[this->curr_index()]; }

      void
      clear()
      {
	this->m_curr_index = 0;
	this->m_max_depth = 0;
	this->m_ival.clear();
      }

      /**
       * Push a new segment into the current heap.
       * N.B. this->begin() includes curr_index()!
       */
      void
      push(const interval& iv)
      {
	this->m_ival.push_back(iv);
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
	this->m_ival.pop_back();
      }

      /**
       * Return the lower limit for the segment at start + ii.
       */
      Tp
      lower_lim(std::size_t ii = 0) const
      { return this->m_ival[this->curr_index() + ii].lower_lim; }

      /**
       * Return the upper limit for the segment at start + ii.
       */
      Tp
      upper_lim(std::size_t ii = 0) const
      { return this->m_ival[this->curr_index() + ii].upper_lim; }

      /**
       * Return the integration result for the segment at start + ii.
       */
      AreaTp
      result(std::size_t ii = 0) const
      { return this->m_ival[this->curr_index() + ii].result; }

      /**
       * Return the absolute error for the segment at start + ii.
       */
      ErrorTp
      abs_error(std::size_t ii = 0) const
      { return this->m_ival[this->curr_index() + ii].abs_error; }

      /**
       * Return the subdivision depth for the segment at start + ii.
       */
      size_t
      depth(std::size_t ii = 0) const
      { return this->m_ival[this->curr_index() + ii].depth; }

      /**
       * Set the absolute error of the segment at start + ii to the given value.
       * Only used by qagp.
       */
      ErrorTp
      set_abs_error(std::size_t ii, ErrorTp abserr)
      { return this->m_ival[this->curr_index() + ii].abs_error = abserr; }

      /**
       * Set the depth of segment at start + ii to the given value.
       * Only used by qagp.
       */
      void
      set_depth(std::size_t ii, std::size_t d)
      { this->m_ival[this->curr_index() + ii].depth = d; }

      bool increment_curr_index();

      /**
       * Reset the index of the current segment to zero.
       * The entire segment vector will me made into a heap.
       */
      void
      reset_curr_index()
      {
	this->m_curr_index = 0;
	this->sort_error();
      }

      /**
       * Return the depth of the current segment.
       */
      std::size_t
      curr_depth() const
      { return this->m_ival[this->curr_index()].depth; }

      std::size_t
      max_depth() const
      { return this->m_max_depth; }

      /**
       * Return the index of the current segment in the segment vector.
       * The segments in [start, end) will be maintained in a heap.
       */
      std::size_t
      curr_index() const
      { return this->m_curr_index; }

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
      AreaTp
      total_integral() const
      {
	auto result_sum = AreaTp{0};
	for (const auto& iv : this->m_ival)
	  result_sum += iv.result;
	return result_sum;
      }

      /// Return the sum of the absolute errors over all integration segments.
      ErrorTp
      total_error() const
      {
	auto tot_error = ErrorTp{0};
	for (auto& iv : this->m_ival)
	  tot_error += iv.abs_error;
	return tot_error;
      }

      /// Return the vector of integration intervals.
      const std::vector<interval>&
      intervals() const
      { return this->m_ival; }

      /*
       * 
       */
      static bool
      subinterval_too_small(Tp a1, Tp a2, Tp b2)
      {
	const auto s_eps = Tp{100} * std::numeric_limits<Tp>::epsilon();
	const auto s_min = Tp{1000} * std::numeric_limits<Tp>::min();

	const auto tmp = (Tp{1} + s_eps) * (std::abs(a2) + s_min);

	return std::abs(a1) <= tmp
	    && std::abs(b2) <= tmp;
      }
    };

  template<typename Tp, typename RetTp>
    std::ostream&
    operator<<(std::ostream& out,
	       const integration_workspace<Tp, RetTp>& ws);

} // namespace emsr

#include <ext/integration_workspace.tcc>

#endif // INTEGRATION_WORKSPACE_H

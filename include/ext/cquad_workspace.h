/* integration/cquad_workspace.h
 *
 * Copyright (C) 2010 Pedro Gonnet
 * Copyright (C) 2016-2020 Free Software Foundation, Inc.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 3 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
 */
// Ported from GSL by Ed Smith-Rowland
// Originally written by Pedro Gonnet
//
// This file contains constants for use in integration schemes.
// Based on structs in gsl/integration/gsl_integration.h

#ifndef CQUAD_WORKSPACE_H
#define CQUAD_WORKSPACE_H 1

#include <vector>

namespace __gnu_cxx
{

  /**
   * Data of a single interval.
   */
  template<typename Tp, typename RetTp>
    struct cquad_interval
    {
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      Tp m_lower_lim;
      Tp m_upper_lim;
      AreaTp m_result;
      AbsAreaTp m_abs_error;
      RetTp m_coeff[64];
      std::array<RetTp, 33> fx;
      std::size_t depth;
      std::size_t rdepth;
      std::size_t ndiv;
    };

  /**
   * Comparison of cquad intervals.
   */
  template<typename Tp, typename RetTp>
    bool
    operator<(const cquad_interval<Tp, RetTp>& ivl,
	      const cquad_interval<Tp, RetTp>& ivr)
    { return ivl.m_abs_error < ivr.m_abs_error; }

  template<typename Tp, typename RetTp>
    struct cquad_interval_comp
    {
      bool
      operator()(const cquad_interval<Tp, RetTp>& ivl,
		 const cquad_interval<Tp, RetTp>& ivr)
      { return ivl.m_abs_error < ivr.m_abs_error; }
    };

  /**
   * The workspace is a collection of intervals.
   * Actually, it is a priority queue where the priority
   * is the absolute error of the interval integral.
   */
  template<typename Tp, typename RetTp>
    struct cquad_workspace
    {
      using AreaTp = decltype(RetTp{} * Tp{});
      using AbsAreaTp = decltype(std::abs(AreaTp{}));

      std::vector<cquad_interval<Tp, RetTp>> m_ival;

      cquad_workspace(std::size_t len = 200)
      : m_ival()
      { this->m_ival.reserve(len); }

      std::size_t size() const
      { return this->m_ival.size(); }

      std::size_t
      capacity() const
      { return this->m_ival.capacity(); }

      typename std::vector<cquad_interval<Tp, RetTp>>::iterator
      begin()
      { return this->m_ival.begin(); }

      typename std::vector<cquad_interval<Tp, RetTp>>::iterator
      end()
      { return this->m_ival.end(); }

      const cquad_interval<Tp, RetTp>&
      top() const
      { return this->m_ival[0]; }

      cquad_interval<Tp, RetTp>&
      top()
      { return this->m_ival[0]; }

      void
      clear()
      { this->m_ival.clear(); }

      void
      push(const cquad_interval<Tp, RetTp>& iv)
      {
	this->m_ival.push_back(iv);
	std::push_heap(this->begin(), this->end(),
		       cquad_interval_comp<Tp, RetTp>{});
      }

      void
      pop()
      {
	std::pop_heap(this->begin(), this->end());
	this->m_ival.pop_back();
      }

      void
      update()
      {
	std::make_heap(this->begin(), this->end(),
		       cquad_interval_comp<Tp, RetTp>{});
      }

      AreaTp
      total_integral() const
      {
	auto tot_igral = AreaTp{0};
	for (auto& iv : m_ival)
	  tot_igral += iv.m_result;
	return tot_igral;
      }

      AbsAreaTp
      total_error() const
      {
	auto tot_error = AbsAreaTp{0};
	for (auto& iv : m_ival)
	  tot_error += iv.m_abs_error;
	return tot_error;
      }
    };

} // namespace __gnu_cxx

#endif // CQUAD_WORKSPACE_H

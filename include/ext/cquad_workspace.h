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
  template<typename _Tp, typename _RetTp>
    struct cquad_interval
    {
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _AbsAreaTp = decltype(std::abs(_AreaTp{}));

      _Tp _M_lower_lim;
      _Tp _M_upper_lim;
      _AreaTp _M_result;
      _AbsAreaTp _M_abs_error;
      _RetTp _M_coeff[64];
      std::array<_RetTp, 33> fx;
      std::size_t depth;
      std::size_t rdepth;
      std::size_t ndiv;
    };

  /**
   * Comparison of cquad intervals.
   */
  template<typename _Tp, typename _RetTp>
    bool
    operator<(const cquad_interval<_Tp, _RetTp>& __ivl,
	      const cquad_interval<_Tp, _RetTp>& __ivr)
    { return __ivl._M_abs_error < __ivr._M_abs_error; }

  template<typename _Tp, typename _RetTp>
    struct cquad_interval_comp
    {
      bool
      operator()(const cquad_interval<_Tp, _RetTp>& __ivl,
		 const cquad_interval<_Tp, _RetTp>& __ivr)
      { return __ivl._M_abs_error < __ivr._M_abs_error; }
    };

  /**
   * The workspace is a collection of intervals.
   * Actually, it is a priority queue where the priority
   * is the absolute error of the interval integral.
   */
  template<typename _Tp, typename _RetTp>
    struct cquad_workspace
    {
      using _AreaTp = decltype(_RetTp{} * _Tp{});
      using _AbsAreaTp = decltype(std::abs(_AreaTp{}));

      std::vector<cquad_interval<_Tp, _RetTp>> _M_ival;

      cquad_workspace(std::size_t __len = 200)
      : _M_ival()
      { this->_M_ival.reserve(__len); }

      std::size_t size() const
      { return this->_M_ival.size(); }

      std::size_t
      capacity() const
      { return this->_M_ival.capacity(); }

      typename std::vector<cquad_interval<_Tp, _RetTp>>::iterator
      begin()
      { return this->_M_ival.begin(); }

      typename std::vector<cquad_interval<_Tp, _RetTp>>::iterator
      end()
      { return this->_M_ival.end(); }

      const cquad_interval<_Tp, _RetTp>&
      top() const
      { return this->_M_ival[0]; }

      cquad_interval<_Tp, _RetTp>&
      top()
      { return this->_M_ival[0]; }

      void
      clear()
      { this->_M_ival.clear(); }

      void
      push(const cquad_interval<_Tp, _RetTp>& __iv)
      {
	this->_M_ival.push_back(__iv);
	std::push_heap(this->begin(), this->end(),
		       cquad_interval_comp<_Tp, _RetTp>{});
      }

      void
      pop()
      {
	std::pop_heap(this->begin(), this->end());
	this->_M_ival.pop_back();
      }

      void
      update()
      {
	std::make_heap(this->begin(), this->end(),
		       cquad_interval_comp<_Tp, _RetTp>{});
      }

      _AreaTp
      total_integral() const
      {
	auto __tot_igral = _AreaTp{0};
	for (auto& __iv : _M_ival)
	  __tot_igral += __iv._M_result;
	return __tot_igral;
      }

      _AbsAreaTp
      total_error() const
      {
	auto __tot_error = _AbsAreaTp{0};
	for (auto& __iv : _M_ival)
	  __tot_error += __iv._M_abs_error;
	return __tot_error;
      }
    };

} // namespace __gnu_cxx

#endif // CQUAD_WORKSPACE_H

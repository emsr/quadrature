// -*- C++ -*-
// Integration utilities for C++.
//
// Copyright (C) 2011-2018 Free Software Foundation, Inc.
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

#ifndef INTEGRATION_WORKSPACE_TCC
#define INTEGRATION_WORKSPACE_TCC 1

#include <iostream>
#include <iomanip>

#include "integration_error.h"

namespace __gnu_cxx
{

  /**
   * Rebuild the current heap.
   * N.B. this->begin() includes curr_index()!
   */
  template<typename _Tp, typename _RetTp>
    void
    integration_workspace<_Tp, _RetTp>::sort_error()
    {
      std::make_heap(this->begin(), this->end(), interval_comp{});
      return;
    }

  /**
   *
   */
  template<typename _Tp, typename _RetTp>
    void
    integration_workspace<_Tp, _RetTp>::
    append(_Tp __a, _Tp __b,
	   _AreaTp __area, _ErrorTp __error,
	   std::size_t __depth)
    {
      interval __iv;
      __iv.__lower_lim = __a;
      __iv.__upper_lim = __b;
      __iv.__result = __area;
      __iv.__abs_error = __error;
      __iv.__depth = __depth;
      this->push(__iv);
    }

  /**
   *
   */
  template<typename _Tp, typename _RetTp>
    void
    integration_workspace<_Tp, _RetTp>::
    split(_Tp __ab,
	  _AreaTp __area1, _ErrorTp __error1,
	  _AreaTp __area2, _ErrorTp __error2)
    {
      auto __iv = this->top();
      const auto __a1 = __iv.__lower_lim;
      const auto __b1 = __ab;
      const auto __a2 = __ab;
      const auto __b2 = __iv.__upper_lim;
      const auto __depth = __iv.__depth + 1;
      this->pop();

      interval __iv1;
      __iv1.__lower_lim = __a1;
      __iv1.__upper_lim = __b1;
      __iv1.__result = __area1;
      __iv1.__abs_error = __error1;
      __iv1.__depth = __depth;
      this->push(__iv1);

      interval __iv2;
      __iv2.__lower_lim = __a2;
      __iv2.__upper_lim = __b2;
      __iv2.__result = __area2;
      __iv2.__abs_error = __error2;
      __iv2.__depth = __depth;
      this->push(__iv2);

      if (__depth > this->_M_max_depth)
	this->_M_max_depth = __depth;
    }

  /**
   * Increase the heap start point until the current segment has a smaller
   * depth than the current maximum depth.  After each increment rebuild
   * the heap from the new start point so the new start point is the largest
   * error.
   *
   * Usage:
   * In the caller the smallest interval (at max depth) has the largest error.
   * Before bisecting decrease the sum of the errors over the larger intervals
   * (error_over_large_intervals) and perform extrapolation.
   */
  template<typename _Tp, typename _RetTp>
    bool
    integration_workspace<_Tp, _RetTp>::increment_curr_index()
    {
      size_t __limit = this->max_size();
      size_t __last = this->size() - 1 ;
      size_t __jupbnd = __last > 1 + __limit / 2
		      ? __limit + 1 - __last
		      : __last;

      const auto __i_max = this->curr_index();
      for (auto __k = __i_max; __k <= __jupbnd; ++__k)
	{
	  if (this->_M_ival[__k].__depth < this->max_depth())
	    return true;
	  else if (this->curr_index() + 1 < this->size())
	    {
	      ++this->_M_curr_index;
	      this->sort_error();
	    }
	}
      return false;
    }

  /**
   * Output the integration workspace to a stream.
   */
  template<typename _Tp, typename _RetTp>
    std::ostream&
    operator<<(std::ostream& __out,
	       const integration_workspace<_Tp, _RetTp>& __ws)
    {
      auto __w = __out.width();
      __out << std::setw(0);
      __out << ' ' << std::setw(2) << __ws.max_depth() << '\n';
      __out << ' ' << std::setw(2) << __ws.curr_index() << '\n';
      for (const auto& __seg : __ws.intervals())
	__out << ' ' << std::setw(2) << __seg.__depth
	      << ' ' << std::setw(__w) << __seg.__lower_lim
	      << ' ' << std::setw(__w) << __seg.__upper_lim
	      << ' ' << std::setw(__w) << __seg.__result
	      << ' ' << std::setw(__w) << __seg.__abs_error
	      << '\n';
      return __out;
    }

} // namespace __gnu_cxx

#endif // INTEGRATION_WORKSPACE_TCC

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

#ifndef INTEGRATION_WORKSPACE_TCC
#define INTEGRATION_WORKSPACE_TCC 1

#include <iostream>
#include <iomanip>

#include <emsr/integration_error.h>

namespace emsr
{

  /**
   * Rebuild the current heap.
   * N.B. this->begin() includes curr_index()!
   */
  template<typename Tp, typename RetTp>
    void
    integration_workspace<Tp, RetTp>::sort_error()
    {
      std::make_heap(this->begin(), this->end(), interval_comp{});
      return;
    }

  /**
   *
   */
  template<typename Tp, typename RetTp>
    void
    integration_workspace<Tp, RetTp>::
    append(Tp a, Tp b,
	   AreaTp area, ErrorTp error,
	   std::size_t depth)
    {
      interval iv;
      iv.lower_lim = a;
      iv.upper_lim = b;
      iv.result = area;
      iv.abs_error = error;
      iv.depth = depth;
      this->push(iv);
    }

  /**
   *
   */
  template<typename Tp, typename RetTp>
    void
    integration_workspace<Tp, RetTp>::
    split(Tp ab,
	  AreaTp area1, ErrorTp error1,
	  AreaTp area2, ErrorTp error2)
    {
      auto iv = this->top();
      const auto a1 = iv.lower_lim;
      const auto b1 = ab;
      const auto a2 = ab;
      const auto b2 = iv.upper_lim;
      const auto depth = iv.depth + 1;
      this->pop();

      interval iv1;
      iv1.lower_lim = a1;
      iv1.upper_lim = b1;
      iv1.result = area1;
      iv1.abs_error = error1;
      iv1.depth = depth;
      this->push(iv1);

      interval iv2;
      iv2.lower_lim = a2;
      iv2.upper_lim = b2;
      iv2.result = area2;
      iv2.abs_error = error2;
      iv2.depth = depth;
      this->push(iv2);

      if (depth > this->m_max_depth)
	this->m_max_depth = depth;
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
  template<typename Tp, typename RetTp>
    bool
    integration_workspace<Tp, RetTp>::increment_curr_index()
    {
      size_t limit = this->max_size();
      size_t last = this->size() - 1 ;
      size_t jupbnd = last > 1 + limit / 2
		      ? limit + 1 - last
		      : last;

      const auto i_max = this->curr_index();
      for (auto k = i_max; k <= jupbnd; ++k)
	{
	  if (this->m_ival[k].depth < this->max_depth())
	    return true;
	  else if (this->curr_index() + 1 < this->size())
	    {
	      ++this->m_curr_index;
	      this->sort_error();
	    }
	}
      return false;
    }

  /**
   * Output the integration workspace to a stream.
   */
  template<typename Tp, typename RetTp>
    std::ostream&
    operator<<(std::ostream& out,
	       const integration_workspace<Tp, RetTp>& ws)
    {
      auto w = out.width();
      out << std::setw(0);
      out << ' ' << std::setw(2) << ws.max_depth() << '\n';
      out << ' ' << std::setw(2) << ws.curr_index() << '\n';
      for (const auto& seg : ws.intervals())
	out << ' ' << std::setw(2) << seg.depth
	      << ' ' << std::setw(w) << seg.lower_lim
	      << ' ' << std::setw(w) << seg.upper_lim
	      << ' ' << std::setw(w) << seg.result
	      << ' ' << std::setw(w) << seg.abs_error
	      << '\n';
      return out;
    }

} // namespace emsr

#endif // INTEGRATION_WORKSPACE_TCC

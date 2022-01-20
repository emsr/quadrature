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
// Ported from GSL by Jason Dick and Ed Smith-Rowland
// Originally written by Brian Gaugh
//
// This file implements an extrapolation table for use in integration schemes
// Based on gsl/integration/qelg.c

#ifndef EXTRAPOLATION_TABLE_H
#define EXTRAPOLATION_TABLE_H 1

#include <array>
#include <utility>
#include <limits>
#include <cmath>
#include <algorithm>

namespace emsr
{

  template<typename AreaTp, typename AbsAreaTp>
    class extrapolation_table
    {

    private:

      std::size_t m_nn;
      std::array<AreaTp, 52> m_rlist2;
      std::size_t m_nres;
      std::array<AreaTp, 3> m_res3la;

      const AbsAreaTp m_irreg_test = AbsAreaTp{0.0001};

    public:

      extrapolation_table()
      : m_nn(0),
	m_nres(0)
      {
	// Try to adjust tests for varing precision?
	// I don't think this is a precision-sensitive limit.
	//this->m_irreg_test = std::pow(AbsAreaTp{10},
	//			 -std::numeric_limits<AbsAreaTp>::digits10 / AbsAreaTp{4.0});
      }

      explicit extrapolation_table(AreaTp y)
      : m_nn(0),
	m_nres(0)
      { this->append(y); }

      void
      append(AreaTp y)
      {
	if (this->m_nn < this->m_rlist2.size())
	  {
	    this->m_rlist2[this->m_nn] = y;
	    ++this->m_nn;
	  }
      }

      std::tuple<AreaTp, AbsAreaTp> qelg();

      std::size_t
      get_nn() const
      { return this->m_nn; }
    };

} // namespace emsr

#include <emsr/extrapolation_table.tcc>

#endif // EXTRAPOLATION_TABLE_H

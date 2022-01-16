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

namespace __gnu_cxx
{

  template<typename _AreaTp, typename _AbsAreaTp>
    class extrapolation_table
    {

    private:

      std::size_t _M_nn;
      std::array<_AreaTp, 52> _M_rlist2;
      std::size_t _M_nres;
      std::array<_AreaTp, 3> _M_res3la;

      const _AbsAreaTp _M_irreg_test = _AbsAreaTp{0.0001};

    public:

      extrapolation_table()
      : _M_nn(0),
	_M_nres(0)
      {
	// Try to adjust tests for varing precision?
	// I don't think this is a precision-sensitive limit.
	//this->_M_irreg_test = std::pow(_AbsAreaTp{10},
	//			 -std::numeric_limits<_AbsAreaTp>::digits10 / _AbsAreaTp{4.0});
      }

      explicit extrapolation_table(_AreaTp __y)
      : _M_nn(0),
	_M_nres(0)
      { this->append(__y); }

      void
      append(_AreaTp __y)
      {
	if (this->_M_nn < this->_M_rlist2.size())
	  {
	    this->_M_rlist2[this->_M_nn] = __y;
	    ++this->_M_nn;
	  }
      }

      std::tuple<_AreaTp, _AbsAreaTp> qelg();

      std::size_t
      get_nn() const
      { return this->_M_nn; }
    };

} // namespace __gnu_cxx

#include <ext/extrapolation_table.tcc>

#endif // EXTRAPOLATION_TABLE_H

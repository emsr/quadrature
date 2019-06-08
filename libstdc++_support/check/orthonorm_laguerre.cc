// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
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

#include <limits>
#include <cassert>
#include <cmath>
#include "integration.h"

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_laguerre(int n1, int n2, _Tp x)
  {
    return std::exp(-x)
	 * std::laguerre(n1, x)
	 * std::laguerre(n2, x);
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_assoc_laguerre()
  {
    const _Tp eps = std::numeric_limits<_Tp>::epsilon();
    const _Tp integ_prec = _Tp{100} * eps;
    const _Tp cmp_prec = _Tp{10} * integ_prec;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1 = itop, n2](_Tp x)
			-> _Tp
			{ return norm_laguerre<_Tp>(n1, n2, x); };

	    auto [result, error]
		= __gnu_cxx::integrate_to_infinity(func, _Tp{0},
						   integ_prec, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_assoc_laguerre<float>();

  test_assoc_laguerre<double>();

  test_assoc_laguerre<long double>();
}

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
#include <ext/cmath>
#include "integration.h"

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename _Tp>
  _Tp
  norm_chebyshev_u(int n1, int n2, _Tp x)
  {
    const auto _S_pi_2 = __gnu_cxx::__math_constants<_Tp>::__pi_half;
    return __gnu_cxx::chebyshev_u(n2, x)
	 * __gnu_cxx::chebyshev_u(n1, x)
	 * std::sqrt(_Tp{1} - x * x)
	 / _S_pi_2;
  }

template<typename _Tp>
  _Tp
  delta(int n1, int n2)
  { return n1 == n2 ? _Tp{1} : _Tp{0}; }

template<typename _Tp>
  void
  test_chebyshev_u()
  {
    const auto eps = std::numeric_limits<_Tp>::epsilon();
    const auto integ_prec = _Tp{100} * eps;
    const auto cmp_prec = _Tp{10} * integ_prec;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](_Tp x)
			-> _Tp
			{ return norm_chebyshev_u(n1, n2, x); };

	    auto [result, error]
		= __gnu_cxx::integrate(func, _Tp{-1}, _Tp{1},
				       integ_prec, _Tp{0});

	    assert(std::abs(delta<_Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_chebyshev_u<float>();

  test_chebyshev_u<double>();

  test_chebyshev_u<long double>();
}

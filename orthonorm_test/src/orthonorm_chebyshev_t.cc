// -*- C++ -*-
// Integration utilities for the C++ library testsuite.
//
// Copyright (C) 2016-2019 Free Software Foundation, Inc.
// Copyright (C) 2020-2022 Edward M. Smith-Rowland
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

#include <limits>
#include <cassert>
#include <cmath>

#include <ext/integration.h>

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_chebyshev_t(int n1, int n2, Tp x)
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto inf = std::numeric_limits<Tp>::infinity();
    if (std::abs(x - Tp{1}) < eps)
      return inf;
    else if (std::abs(x + Tp{1}) < eps)
      return (n1 + n2) & 1 ? -inf : inf;
    else
      return __gnu_cxx::chebyshev_t(n2, x)
	   * __gnu_cxx::chebyshev_t(n1, x)
	   / std::sqrt(Tp{1} - x * x);
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_chebyshev_t()
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto integ_prec = Tp{100} * eps;
    constexpr auto cmp_prec = Tp{10} * integ_prec;

    for (int n1 : {0, 5, 10})
      {
	for (int n2 : {0, 5, 10})
	  {
	    auto func = [n1, n2](Tp x)
			-> Tp
			{ return norm_chebyshev_t(n1, n2, x); };

	    auto [result, error]
		= emsr::integrate_singular_endpoints(func,
				 Tp{-1}, Tp{1},
				 Tp{-0.5}, Tp{-0.5}, 0, 0,
				 integ_prec, Tp{0});

	    assert(std::abs(delta<Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_chebyshev_t<float>();

  test_chebyshev_t<double>();

  test_chebyshev_t<long double>();
}

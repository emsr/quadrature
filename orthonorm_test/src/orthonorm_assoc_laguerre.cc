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

// Try to manage the gamma ratio.
template<typename Tp>
  Tp
  gamma_ratio(int n, int m)
  {
    auto gaman1 = std::tgamma(Tp(1 + m));
    auto fact = gaman1;
    for (int k = 1; k <= n; ++k)
      fact *= Tp(k + m) / Tp(k);
    return fact;
  }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  norm_assoc_laguerre(int n1, int n2, int m, Tp x)
  {
    auto norm = gamma_ratio<Tp>(n1, m);
    return std::pow(x, m) * std::exp(-x)
	 * std::assoc_laguerre(n1, m, x)
	 * std::assoc_laguerre(n2, m, x) / norm;
  }

template<typename Tp>
  Tp
  delta(int n1, int n2)
  { return n1 == n2 ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_assoc_laguerre(Tp alpha)
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto integ_prec = Tp{1000} * eps;
    constexpr auto cmp_prec = Tp{10} * integ_prec;

    for (int n1 : {0, 2, 5})
      {
	for (int n2 : {0, 2, 5})
	  {
	    auto func = [n1, n2, alpha](Tp x)
			-> Tp
			{ return norm_assoc_laguerre<Tp>(n1, n2, alpha, x); };

	    auto [result, error]
		= emsr::integrate_lower_pinf(func, Tp{0}, integ_prec, Tp{0});

	    assert(std::abs(delta<Tp>(n1, n2) - result) < cmp_prec);
	  }
      }
  }

int
main()
{
  test_assoc_laguerre<float>(0);

  test_assoc_laguerre<double>(0);

  test_assoc_laguerre<long double>(0);
}

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

#include <emsr/integration.h>

/**
 * Neumann's number
 */
template<typename Tp>
  Tp
  epsilon(int m)
  { return m == 0 ? Tp{2} : Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
template<typename Tp>
  Tp
  normalized_radpoly(int n1, int m1, int n2, int m2, Tp rho)
  {
    auto norm = Tp{1} / std::sqrt(Tp(2 * n1 + 2) * Tp(2 * n2 + 2));
    return rho
	 * emsr::radpoly(n1, m1, rho)
	 * emsr::radpoly(n2, m2, rho) / norm;
  }

template<typename Tp>
  Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_radpoly()
  {
    constexpr auto eps = std::numeric_limits<Tp>::epsilon();
    constexpr auto integ_precision = Tp{1000} * eps;
    constexpr auto comp_precision = Tp{10} * integ_precision;

    for (int n1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    if ((n1 - m1) & 1)
	      continue;
	    for (int n2 : {0, 2, 5})
	      {
		// The orthonormality only works for m2 == m1.
		int m2 = m1;
		if (m2 > n2)
		  continue;
		//for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;

		    auto func = [n1, m1, n2, m2](Tp x)
				-> Tp
			      { return normalized_radpoly(n1, m1, n2, m2, x); };

		    auto [result, error]
			= emsr::integrate_singular(func, Tp{0}, Tp{1},
					     integ_precision, Tp{0});

		    assert(std::abs(delta<Tp>(n1, n2) - result) < cmp_prec);
		  }
	      }
	  }
      }
  }

int
main()
{
  test_radpoly<float>();

  test_radpoly<double>();

  test_radpoly<long double>();
}

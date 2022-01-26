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
#include <numbers>

#include <emsr/integration.h>
#include <emsr/math_constants.h>

/**
 * Neumann's number
 */
template<typename Tp>
  Tp
  epsilon(int m)
  { return m == 0 ? Tp{2} : Tp{1}; }

// Function which should integrate to 1 for n1 == n2, 0 otherwise.
// This actually does the angular integral.
template<typename Tp>
  Tp
  norm_zernike(int n1, int m1, int n2, int m2, Tp rho)
  {
    constexpr auto eps = Tp{10000} * std::numeric_limits<Tp>::epsilon();
    constexpr auto 2pi = Tp{2} * emsr::pi_v<Tp>;
    auto z1 = [n1, m1, rho](Tp phi)
	      ->Tp { return emsr::zernike(n1, m1, rho, phi); };
    auto z2 = [n2, m2, rho](Tp phi)
	      ->Tp { return emsr::zernike(n2, m2, rho, phi); };
    auto norm = Tp{1} / std::sqrt(Tp(2 * n1 + 2) * Tp(2 * n2 + 2));
    auto fun = [n1, m1, rho, z1, z2, norm](Tp phi)
	       ->Tp { return rho * z1(phi) * z2(phi) / norm; };
    Tp val;
    std::tie(val, std::ignore) = emsr::integrate(fun, Tp{0}, Tp{2pi},
					   eps, Tp{0}, 1024, QK_61);
    return Tp{2} * val / 2pi / epsilon<Tp>(m1);
  }

template<typename Tp>
  Tp
  delta(int n1, int m1, int n2, int m2)
  { return (n1 == n2 && m1 == m2) ? Tp{1} : Tp{0}; }

template<typename Tp>
  void
  test_zernike()
  {
    const auto eps = std::numeric_limits<Tp>::epsilon();
    const auto integ_prec = Tp{1000} * eps;
    const auto cmp_prec = Tp{10} * integ_prec;

    for (int n1 : {0, 2, 5})
      {
	for (int m1 = 0; m1 <= n1; ++m1)
	  {
	    if ((n1 - m1) & 1)
	      continue;
	    for (int n2 : {0, 2, 5})
	      {
		for (int m2 = 0; m2 <= n2; ++m2)
		  {
		    if ((n2 - m2) & 1)
		      continue;
		    auto func = [n1, m1, n2, m2](Tp x)
				-> Tp
				{ return norm_zernike(n1, m1, n2, m2, x); };

		    auto [result, error]
			= emsr::integrate_singular(func, Tp{0}, Tp{1}, integ_prec, Tp{0});

		    assert(std::abs(delta<Tp>(n1, m1, n2, m2) - result)
				 < cmp_prec);
	      }
	  }
      }
  }

int
main()
{
  test_zernike<float>();

  test_zernike<double>();

  test_zernike<long double>();
}
